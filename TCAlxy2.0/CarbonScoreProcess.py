import threading
import numpy as np
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import os

# 线程本地存储用于数据库连接
thread_local = threading.local()

class GreedyMatcher:
    def __init__(self, db_path, user_data, c_range=(10,30), fw_range=(100.0,300.0), mode='typed'):
        """
        :param mode: 匹配模式 ('typed'按类型匹配/'untyped'全局匹配)
        """
        self.db_path = db_path
        self.user_data = user_data
        self.c_range = c_range
        self.fw_range = fw_range
        self.mode = mode
        self._validate_inputs()
        self._init_results_table()
        self._create_indexes()
        self.carbon_cache = {}

    def _validate_inputs(self):
        """参数校验（支持两种模式）"""
        if self.mode == 'typed':
            if not isinstance(self.user_data, dict) or any(not isinstance(v, list) for v in self.user_data.values()):
                raise ValueError("The type pattern requires {0:[],1:[],...} format")
            self.user_total = sum(len(v) for v in self.user_data.values())
        elif self.mode == 'untyped':
            if not isinstance(self.user_data, list):
                raise ValueError("The global mode requires a list of displacements")
            self.user_total = len(self.user_data)
        else:
            raise ValueError("Invalid mode: typed/untyped")

        if self.c_range[0] > self.c_range[1] or self.fw_range[0] > self.fw_range[1]:
            raise ValueError("The filter range is invalid")
        if self.user_total <= 0:
            raise ValueError("User data cannot be empty")

    @contextmanager
    def _get_conn(self):
        """线程安全的数据库连接管理"""
        if not hasattr(thread_local, "conn"):
            thread_local.conn = sqlite3.connect(
                self.db_path,
                check_same_thread=False,
                timeout=30,
                isolation_level=None
            )
            thread_local.conn.execute("PRAGMA journal_mode=WAL")
            thread_local.conn.execute("PRAGMA synchronous=NORMAL")
            thread_local.conn.execute("PRAGMA cache_size=-20000")  # 20MB缓存
        try:
            yield thread_local.conn
        except sqlite3.DatabaseError as e:
            thread_local.conn.rollback()
            raise

    def _create_indexes(self):
        """创建必要索引"""
        with self._get_conn() as conn:
            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_molecules_filter 
                ON molecules(carbon_count, fw)
            """)
            conn.execute("""
                CREATE INDEX IF NOT EXISTS idx_carbons_mol 
                ON carbons(mol_id)
            """)

    def _init_results_table(self):
        """初始化结果表"""
        with self._get_conn() as conn:
            conn.execute("DROP TABLE IF EXISTS greedy_matches")
            conn.execute("""
                CREATE TABLE IF NOT EXISTS greedy_matches (
                    match_id INTEGER PRIMARY KEY,
                    mol_id INTEGER,
                    mol_name TEXT,
                    user_shift REAL,
                    db_c_index INTEGER,
                    db_shift REAL,
                    type_matched BOOLEAN,
                    connected_h_shifts TEXT,
                    connected_carbons TEXT,
                    FOREIGN KEY(mol_id) REFERENCES molecules(mol_id)
                )
            """)

    def _precache_molecules(self):
        """预加载分子数据到内存（修正过滤方向）"""
        with self._get_conn() as conn:
            query = """
                SELECT m.mol_id, m.name, c.c_index, c.h_count, c.c_shift,
                       c.h_shifts, c.connected_carbon
                FROM molecules m
                JOIN carbons c ON m.mol_id = c.mol_id
                WHERE m.carbon_count BETWEEN ? AND ?    -- 用户指定范围
                  AND m.fw BETWEEN ? AND ?              -- 分子量范围
                  AND m.carbon_count <= ?               -- 关键修正点：数据库碳数≤用户总碳数
            """
            params = (
                self.c_range[0],
                self.c_range[1],
                self.fw_range[0],
                self.fw_range[1],
                self.user_total  # 用户总碳数
            )
            rows = conn.execute(query, params).fetchall()

            for row in rows:
                mol_id, name, c_idx, hc, c_shift, h_shifts, conn_c = row
                if mol_id not in self.carbon_cache:
                    self.carbon_cache[mol_id] = {
                        'name': name,
                        'carbons': {}
                    }
                self.carbon_cache[mol_id]['carbons'][c_idx] = {
                    'c_index': c_idx,  # 新增此行
                    'h_count': hc,
                    'c_shift': c_shift,
                    'h_shifts': h_shifts,
                    'connected_c': conn_c
                }

    # def _vectorized_match(self, user_shifts, db_candidates):
    #     """数据库原子优先的向量化匹配，确保每个原子都有匹配"""
    #     if not user_shifts or not db_candidates:
    #         return [], db_candidates
    # 
    #     try:
    #         # 构建差异矩阵（数据库原子为行，用户位移为列）
    #         db = np.array([c['c_shift'] for c in db_candidates], dtype=np.float32)
    #         user = np.array(user_shifts, dtype=np.float32)
    #         diff = np.abs(db[:, np.newaxis] - user)
    # 
    #         matches = []
    #         used_user = set()
    # 
    #         # 为每个数据库原子寻找最佳可用用户位移
    #         for db_idx in range(len(db_candidates)):
    #             if diff[db_idx].min() == np.inf:
    #                 continue
    # 
    #             # 在未使用的用户位移中找最小值
    #             valid_user = [i for i in range(len(user)) if i not in used_user]
    #             if not valid_user:
    #                 break
    # 
    #             u_idx = valid_user[np.argmin(diff[db_idx][valid_user])]
    #             matches.append((
    #                 user_shifts[u_idx],
    #                 db_candidates[db_idx]['c_index'],  # 数据库原子索引
    #                 db_candidates[db_idx]['c_shift']
    #             ))
    #             used_user.add(u_idx)
    #             diff[:, u_idx] = np.inf  # 标记该用户位移为已用
    # 
    #         return matches, []
    # 
    #     except Exception as e:
    #         print(f"匹配错误: {str(e)}")
    #         return [], db_candidates
    def _vectorized_match(self, user_shifts, db_candidates):
        """数据库原子优先的向量化匹配，相同c_shift的原子匹配同一个最佳用户位移"""
        if not user_shifts or not db_candidates:
            return [], db_candidates

        try:
            # 按c_shift分组，相同c_shift的原子归为一组
            shift_groups = {}
            for carbon in db_candidates:
                c_shift = carbon['c_shift']
                if c_shift not in shift_groups:
                    shift_groups[c_shift] = []
                shift_groups[c_shift].append(carbon)

            # 构建差异矩阵（按组处理，每组一个代表）
            group_rep = []  # 每组代表: (c_shift, 组内原子列表)
            for c_shift, carbons in shift_groups.items():
                group_rep.append((c_shift, carbons))

            db_shifts = np.array([rep[0] for rep in group_rep], dtype=np.float32)
            user = np.array(user_shifts, dtype=np.float32)
            diff = np.abs(db_shifts[:, np.newaxis] - user)

            matches = []
            used_user = set()

            # 为每组寻找最佳可用用户位移
            for group_idx in range(len(group_rep)):
                if diff[group_idx].min() == np.inf:
                    continue

                # 在未使用的用户位移中找最小值
                valid_user = [i for i in range(len(user)) if i not in used_user]
                if not valid_user:
                    break

                u_idx = valid_user[np.argmin(diff[group_idx][valid_user])]
                c_shift_val = db_shifts[group_idx]
                group_carbons = group_rep[group_idx][1]

                # 整组原子使用同一个匹配的用户位移
                for carbon in group_carbons:
                    matches.append((
                        user_shifts[u_idx],
                        carbon['c_index'],  # 数据库原子索引
                        carbon['c_shift']
                    ))

                used_user.add(u_idx)
                diff[:, u_idx] = np.inf  # 标记该用户位移为已用

            # 找出未匹配的原子（整组未匹配的原子）
            matched_indices = {carbon['c_index'] for _, c_idx, _ in matches for carbon in db_candidates if carbon['c_index'] == c_idx}
            unmatched = [carbon for carbon in db_candidates if carbon['c_index'] not in matched_indices]

            return matches, unmatched

        except Exception as e:
            print(f"匹配错误: {str(e)}")
            return [], db_candidates
    def _process_molecule(self, mol_id):
        """整合两种匹配模式"""
        try:
            if mol_id not in self.carbon_cache:
                return

            mol_data = self.carbon_cache[mol_id]
            all_matches = []

            if self.mode == 'typed':
                # 原有类型匹配逻辑
                used_user_shifts = set()
                for h_type in [3, 2, 1, 0]:
                    user_shifts = self.user_data.get(h_type, [])
                    candidates = [c for c in mol_data['carbons'].values() if c['h_count'] == h_type]

                    matches, _ = self._vectorized_match(
                        [s for s in user_shifts if s not in used_user_shifts],
                        candidates
                    )
                    for u_shift, c_idx, c_shift in matches:
                        all_matches.append( (h_type, u_shift, c_idx, c_shift, True) )
                        used_user_shifts.add(u_shift)

                # 阶段2：剩余原子全局匹配（使用所有未使用的用户位移）
                remaining_candidates = [
                    c for c in mol_data['carbons'].values()
                    if c['c_index'] not in [m[2] for m in all_matches]
                ]
                remaining_user = [
                    s for h_type in [0,1,2,3]
                    for s in self.user_data[h_type]
                    if s not in used_user_shifts
                ]

                if remaining_candidates and remaining_user:
                    matches, _ = self._vectorized_match(remaining_user, remaining_candidates)
                    all_matches.extend(
                        (None, u_shift, c_idx, c_shift, False)
                        for u_shift, c_idx, c_shift in matches
                    )
            elif self.mode == 'untyped':
                # 全局匹配逻辑
                candidates = list(mol_data['carbons'].values())
                matches, _ = self._vectorized_match(self.user_data, candidates)
                for u_shift, c_idx, c_shift in matches:
                    all_matches.append( (None, u_shift, c_idx, c_shift, True) )

            self._bulk_save_matches(mol_id, mol_data['name'], all_matches)


        except Exception as e:
            print(f"Processing molecule {mol_id} failed: {str(e)}")

    def _bulk_save_matches(self, mol_id, mol_name, matches):
        """保存匹配结果到数据库"""
        data = []
        for match in matches:
            h_type, u_shift, c_idx, db_shift, type_matched = match
            carbon = self.carbon_cache[mol_id]['carbons'].get(c_idx, {})
            data.append((
                mol_id,
                mol_name,
                u_shift,
                c_idx,
                db_shift,
                type_matched,
                carbon.get('h_shifts', '[]'),
                carbon.get('connected_c', '[]')
            ))

        if data:
            try:
                with self._get_conn() as conn:
                    conn.executemany("""
                        INSERT INTO greedy_matches 
                        (mol_id, mol_name, user_shift, 
                         db_c_index, db_shift, type_matched, connected_h_shifts, connected_carbons)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """, data)
                    conn.commit()
            except sqlite3.Error as e:
                print(f"Save molecule {mol_id} result failed: {str(e)}")
    def run(self, max_workers=None):
        """执行匹配流程"""
        # 自动设置线程数
        max_workers = max_workers or min(os.cpu_count() * 2, 32)

        print("Pre-loaded molecular data...")
        self._precache_molecules()
        qualified = list(self.carbon_cache.keys())
        print(f"Find {len(qualified)} molecules that are eligible")

        # 分批次处理 (避免内存溢出)
        batch_size = 500
        batches = [qualified[i:i+batch_size] for i in range(0, len(qualified), batch_size)]

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            with tqdm(total=len(qualified), desc="Match progress") as pbar:
                for batch in batches:
                    futures = {
                        executor.submit(self._process_molecule, mol_id): mol_id
                        for mol_id in batch
                    }
                    for future in as_completed(futures):
                        try:
                            future.result()
                        except Exception as e:
                            mid = futures[future]
                            print(f"Molecule {mid} error: {str(e)}")
                        finally:
                            pbar.update(1)

        # 最终验证
        self._validate_results()

    def _validate_results(self):
        """结果验证"""
        with self._get_conn() as conn:
            # 验证过滤条件
            invalid = conn.execute("""
                SELECT COUNT(*) 
                FROM greedy_matches gm
                JOIN molecules m ON gm.mol_id = m.mol_id
                WHERE m.carbon_count > ?
            """, (self.user_total,)).fetchone()[0]

            if invalid > 0:
                raise ValueError(f"{invalid} violation matches found！")

            # 验证数据完整性
            total_matches = conn.execute("SELECT COUNT(*) FROM greedy_matches").fetchone()[0]
            print(f"Verified! A total of {total_matches} valid matching records were generated")


# import sqlite3
# import json
#
# def get_molecule_details(db_path: str, target_mol_id: int = 695):
#     """获取指定分子的碳原子数据"""
#     conn = sqlite3.connect(db_path)
#
#     # 检查分子是否存在
#     cursor = conn.execute("SELECT EXISTS(SELECT 1 FROM molecules WHERE mol_id=?)", (target_mol_id,))
#     if not cursor.fetchone()[0]:
#         print(f"分子ID {target_mol_id} 不存在")
#         return
#
#     # 获取碳原子数据
#     cursor = conn.execute('''
#         SELECT c_index, h_count, h_shifts, c_shift
#         FROM carbons
#         WHERE mol_id = ?
#         ORDER BY c_index
#     ''', (target_mol_id,))
#
#     results = []
#     for row in cursor.fetchall():
#         results.append({
#             "carbon_index": row[0],
#             "h_count": row[1],
#             "h_shifts": json.loads(row[2]),  # 转换为列表
#             "c_shift": row[3]
#         })
#
#     conn.close()
#     return results
#
# # 使用示例
# data = get_molecule_details("chem_data.db")
# if data:
#     for item in data:
#         print(f"碳原子 {item['carbon_index']}:")
#         print(f"  连接H数: {item['h_count']}")
#         print(f"  H位移: {item['h_shifts']}")
#         print(f"  C位移: {item['c_shift']}\n")

import sqlite3
import numpy as np
from typing import List, Tuple, Union

class CarbonScorer:
    def __init__(self, db_path: str):
        """
        :param db_path: 数据库路径
        """
        self.db_path = db_path
        self._validate_db_structure()

    def _validate_db_structure(self):
        """验证数据库表结构"""
        required_tables = ['greedy_matches']
        with sqlite3.connect(self.db_path) as conn:
            tables = conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
            existing_tables = {t[0] for t in tables}
            missing = set(required_tables) - existing_tables
            if missing:
                raise ValueError(f"Required forms are missing：{missing}")

    def create_results_table(self):
        """创建评分结果表"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("DROP TABLE IF EXISTS carbon_scores")
            conn.execute("""
                CREATE TABLE carbon_scores (
                    score_id INTEGER PRIMARY KEY,
                    mol_id INTEGER,
                    mol_name TEXT,
                    c_index INTEGER,
                    db_shift REAL,
                    user_shift REAL,
                    score REAL,
                    color TEXT,
                    connected_carbons TEXT,
                    FOREIGN KEY(mol_id) REFERENCES molecules(mol_id)
                )
            """)
            conn.commit()

    def process_scoring(
            self,
            mode: str,
            global_thresholds: Tuple[float, float] = None,
            fine_ranges: List[Tuple[Tuple[float, float], float, float]] = None
    ):
        """
        执行评分流程

        :param mode: 模式选择 ('global' 或 'fine')
        :param global_thresholds: 全局阈值 (green_threshold, yellow_threshold)
        :param fine_ranges: 精细模式范围列表 [((low, high), green, yellow), ...]
        """
        # 参数校验
        if mode == 'global':
            if not global_thresholds or len(global_thresholds) != 2:
                raise ValueError("The global mode requires two thresholds")
            green, yellow = global_thresholds
            if green >= yellow:
                raise ValueError("The green threshold must be less than the yellow threshold")
        elif mode == 'fine':
            if not fine_ranges:
                raise ValueError("Granular mode requires a list of ranges")
            for r in fine_ranges:
                (low, high), g, y = r
                if low >= high or g >= y:
                    raise ValueError("Invalid scope or threshold")
        else:
            raise ValueError("Invalid mode selection")

        # 处理数据
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute("""
                SELECT gm.mol_id, gm.mol_name, gm.db_c_index as c_index,
                       gm.db_shift, gm.user_shift, gm.type_matched,
                       gm.connected_carbons
                FROM greedy_matches gm
            """)

            batch = []
            for row in cursor:
                # 获取基础数据
                db_shift = row['db_shift']
                user_shift = row['user_shift']
                type_matched = row['type_matched']

                if db_shift is not None and user_shift is not None:
                    delta = abs(db_shift - user_shift) if user_shift else None
                else:
                    delta = 1000

                # 计算评分和颜色
                if mode == 'global':
                    score, color = self._global_score(delta, type_matched, global_thresholds)
                else:
                    score, color = self._fine_score(db_shift, delta, type_matched, fine_ranges)

                batch.append((
                    row['mol_id'],
                    row['mol_name'],
                    row['c_index'],
                    db_shift,
                    user_shift,
                    score,
                    color,
                    row['connected_carbons']
                ))

                # 批量插入
                if len(batch) >= 1000:
                    self._insert_batch(conn, batch)
                    batch = []

            if batch:
                self._insert_batch(conn, batch)

    def _global_score(
            self,
            delta: float,
            type_matched: bool,
            thresholds: Tuple[float, float]
    ) -> Tuple[float, str]:
        """全局模式评分"""
        if not type_matched or delta is None:
            return 0.0, '#FF0000'

        green, yellow = thresholds
        if delta < green:
            return 1.0, '#00FF00'
        elif delta > yellow:
            return 0.0, '#FF0000'
        else:
            # 线性插值
            score = 1 - (delta - green) / (yellow - green)
            return round(score, 2), '#FFFF00'

    def _fine_score(
            self,
            db_shift: float,
            delta: float,
            type_matched: bool,
            ranges: List[Tuple[Tuple[float, float], float, float]]
    ) -> Tuple[float, str]:
        """精细模式评分"""
        if not type_matched or delta is None:
            return 0.0, 'red'

        # 查找匹配范围（按输入顺序优先）
        matched_range = None
        for r in ranges:
            (low, high), g, y = r
            if low <= db_shift <= high:
                matched_range = (g, y)
                break

        if not matched_range:
            return 0.0, 'red'

        green, yellow = matched_range
        if delta < green:
            return 1.0, 'green'
        elif delta > yellow:
            return 0.0, 'red'
        else:
            score = 1 - (delta - green) / (yellow - green)
            return round(score, 2), 'yellow'

    def _insert_batch(self, conn, batch):
        """批量插入数据"""
        conn.executemany("""
            INSERT INTO carbon_scores 
            (mol_id, mol_name, c_index,db_shift,user_shift, score, color, connected_carbons)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, batch)
        conn.commit()

# 使用示例
# if __name__ == "__main__":
#     # 初始化
#     scorer = CarbonScorer("chem_data.db")
#     scorer.create_results_table()
#
#     # # 示例1：全局模式
#     # scorer.process_scoring(
#     #     mode='global',
#     #     global_thresholds=(0.5, 2.0)
#     # )
#
#     # 示例2：精细模式
#     fine_ranges = [
#         ((0, 50), 1, 2.0),
#         ((50, 150), 2, 3.0),
#         ((150, 300), 3, 4.0)
#     ]
#     scorer.process_scoring(
#         mode='fine',
#         fine_ranges=fine_ranges
#     )

import sqlite3
import json
from collections import deque
from typing import Dict, List, Tuple
import logging

# 配置日志记录
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

class RobustMoleculeScorer:
    def __init__(self, db_path: str):
        self.db_path = db_path
        self._validate_structure()

    def _validate_structure(self):
        """验证数据库表结构完整性"""
        required_tables = {'carbon_scores', 'molecules', 'greedy_matches'}
        with sqlite3.connect(self.db_path) as conn:
            existing_tables = {row[0] for row in
                               conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            missing = required_tables - existing_tables
            if missing:
                raise ValueError(f"缺失关键数据表: {missing}")

    def create_result_table(self):
        """创建分子评分结果表"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("DROP TABLE IF EXISTS molecule_scores")
            conn.execute("""
                CREATE TABLE IF NOT EXISTS molecule_scores (
                    mol_id INTEGER PRIMARY KEY,
                    mol_name TEXT,
                    final_score REAL,
                    score_details TEXT,
                    FOREIGN KEY(mol_id) REFERENCES molecules(mol_id)
                )
            """)
            conn.commit()

    def process_molecules(
            self,
            self_weight: float,
            env_weight: float,
            env_level: int = 1
    ):
        """
        增强鲁棒性的评分流程

        :param self_weight: 自身评分权重 (0-1)
        :param env_weight: 环境评分权重 (0-1)
        :param env_level: 环境层级 (1-3)
        """
        # 权重校验
        if not (0 <= self_weight <= 1) or not (0 <= env_weight <= 1):
            raise ValueError("The weight parameter must be between 0 and 1")
        if abs(self_weight + env_weight - 1.0) > 1e-6:
            raise ValueError("The sum of the weights must be equal to 1")
        if env_level not in {1, 2, 3}:
            raise ValueError("The environment level must be 1, 2, or 3")

        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            # 获取所有分子ID（按分子ID排序保证可重复性）
            mol_ids = [row[0] for row in
                       conn.execute("SELECT DISTINCT mol_id FROM greedy_matches")]

            # 初始化进度条
            total_mols = len(mol_ids)
            with tqdm(total=total_mols, desc="Molecular scoring progress", unit="mol") as pbar:
                for mol_id in mol_ids:
                    try:
                        # 加载分子数据（包含严格校验）
                        mol_data = self._load_molecule_data(conn, mol_id)
                        # 构建连接关系图（带有效性检查）
                        connection_graph = self._build_connection_graph(mol_data)
                        # 计算分子评分
                        self._score_molecule(conn, mol_id, mol_data,
                                             connection_graph, self_weight,
                                             env_weight, env_level)
                    except Exception as e:
                        logging.error(f"Processing Molecule {mol_id} Failed: {str(e)}")
                        raise
                    finally:
                        # 无论成功失败都更新进度
                        pbar.update(1)
                        pbar.set_postfix_str(f"Current molecule: {mol_id}")

        print(f"\nGrading done! Processed {total_mols} molecules")


    def _load_molecule_data(self, conn, mol_id: int) -> Dict[int, dict]:
        """加载分子数据并进行完整性校验"""
        cursor = conn.execute("""
            SELECT c_index, score, connected_carbons 
            FROM carbon_scores 
            WHERE mol_id = ?
        """, (mol_id,))

        molecule_data = {}
        for row in cursor:
            # 数据校验
            if not (0 <= row['score'] <= 1):
                raise ValueError(
                    f"The carbon atom {row['c_index']} score of the molecule {mol_id} is out of range: {row['score']}"
                )

            try:
                connections = json.loads(row['connected_carbons'])
                if not isinstance(connections, list):
                    raise TypeError("The connection relationship must be a list")
            except (json.JSONDecodeError, TypeError) as e:
                raise ValueError(f"Carbon atom {row['c_index']} of molecule {mol_id} is malformed") from e

            molecule_data[row['c_index']] = {
                'score': row['score'],
                'connections': connections
            }
        return molecule_data

    def _build_connection_graph(self, mol_data: dict) -> Dict[int, List[int]]:
        """构建带有效性检查的连接关系图"""
        valid_indices = set(mol_data.keys())
        graph = {}

        for c_index, data in mol_data.items():
            valid_connections = []
            for neighbor in data['connections']:
                if neighbor in valid_indices:
                    valid_connections.append(neighbor)
                else:
                    pass
            graph[c_index] = valid_connections
        return graph

    def _get_environment_atoms(
            self,
            root: int,
            graph: Dict[int, List[int]],
            max_level: int
    ) -> List[int]:
        """获取多级环境原子（带连接有效性检查）"""
        visited = set([root])
        queue = deque([(root, 0)])
        environment = []

        while queue:
            current_atom, current_level = queue.popleft()
            if current_level >= max_level:
                continue

            for neighbor in graph.get(current_atom, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    environment.append(neighbor)
                    queue.append((neighbor, current_level + 1))
        return environment

    def _calculate_env_score(
            self,
            env_atoms: List[int],
            mol_data: dict
    ) -> float:
        """计算环境评分（带数据完整性检查）"""
        valid_scores = []
        for atom in env_atoms:
            if atom not in mol_data:
                logging.warning(f" Ignore environmental atoms that do not exist in the molecular data:{atom}")
                continue
            valid_scores.append(mol_data[atom]['score'])

        if not valid_scores:
            return 0.0
        return sum(valid_scores) / len(valid_scores)

    def _score_molecule(
            self,
            conn,
            mol_id: int,
            mol_data: dict,
            graph: dict,
            self_weight: float,
            env_weight: float,
            env_level: int
    ):
        """执行分子评分并存储结果"""
        total_score = 0.0
        score_details = []

        for c_index, data in mol_data.items():
            # 获取环境原子
            env_atoms = self._get_environment_atoms(c_index, graph, env_level)
            # 计算环境评分
            try:
                env_score = self._calculate_env_score(env_atoms, mol_data)
            except ZeroDivisionError:
                env_score = 0.0

            # 计算综合评分
            if data['score'] == 0:
                composite = 0
            else:
                composite = (self_weight * data['score'] +
                             env_weight * env_score)
            total_score += composite

            # 记录评分细节
            score_details.append({
                'c_index': c_index,
                'self_score': data['score'],
                'env_score': round(env_score, 4),
                'env_atoms': env_atoms,
                'composite': round(composite, 4)
            })

        # 计算最终得分
        final_score = total_score / len(mol_data) if mol_data else 0
        # 获取分子名称
        mol_name = conn.execute(
            "SELECT name FROM molecules WHERE mol_id = ?", (mol_id,)
        ).fetchone()[0]

        # 保存结果
        conn.execute("""
            INSERT OR REPLACE INTO molecule_scores 
            (mol_id, mol_name, final_score, score_details)
            VALUES (?, ?, ?, ?)
        """, (
            mol_id,
            mol_name,
            round(final_score, 4),
            json.dumps(score_details)
        ))
        conn.commit()

# 使用示例
# if __name__ == "__main__":
#     # 初始化评分器
#     scorer = RobustMoleculeScorer("chem_data.db")
#     scorer.create_result_table()
#
#     # 执行评分（示例参数）
#     try:
#         scorer.process_molecules(
#             self_weight=0.6,
#             env_weight=0.4,
#             env_level=2
#         )
#     except Exception as e:
#         logging.error(f"评分流程失败: {str(e)}")
#         raise
