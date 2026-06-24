# ChemDraw API — 调用方使用指南（面向 Code Agent）

本服务提供化学名称 / SMILES 互转与分子比较能力。请按本指南调用，以获得稳定、
高效、对服务端友好的使用体验。

---

## 1. 接入信息

- **Base URL**: `https://<服务域名或公网IP>/chemdraw/api`
- **认证**: 除 `/health` 外，所有接口需携带 Bearer Token
  - 请求头: `Authorization: Bearer <API_KEY>`
- **内容类型**: `Content-Type: application/json`
- **字符编码**: UTF-8

> 向服务方索取 `Base URL` 与 `API_KEY`，不要把 key 硬编码进可公开的代码或日志。

---

## 2. 接口一览

| 方法 | 路径 | 用途 | 请求体 | 成功响应 |
|------|------|------|--------|----------|
| GET  | `/health` | 健康检查（无需认证） | — | `{"status":"ok"}` |
| POST | `/name2smiles` | 化学名称 → SMILES | `{"name": "aspirin"}` | `{"smiles": "..."}` |
| POST | `/smiles2name` | SMILES → 化学名称 | `{"smiles": "CCO"}` | `{"name": "..."}` |
| POST | `/mol_compare` | 比较两个分子 | `{"mol1":"...", "mol2":"..."}` | `{"exact_match":bool,"tanimoto":float,"warning":str\|null}` |
| POST | `/batch_mol_compare` | **批量**比较多对分子（首选） | 见 §4 | 见 §4 |

说明：
- `mol_compare` / `batch_mol_compare` 的分子表示由 ChemScript 自动解析，**SMILES、化学名称等均可**。判等基于 InChI，相似度为 Tanimoto。
- `name2smiles` 若无法识别名称，会返回**空字符串** `{"smiles": ""}`（HTTP 200，不报错）。调用方需自行判断空值。
- ⚠️ `/smiles2rdkit` 当前不可用（返回 500），请勿调用。

---

## 3. 健康使用准则（请务必遵守）

1. **能批量就批量**：要比较多对分子时，用 `batch_mol_compare` 一次提交，
   **不要** for-loop 单条调用。这能把 N 次网络往返压成 1 次，是对带宽和延迟最友好的方式。
2. **控制并发**：客户端并发请求数建议从 **20–50** 起步，观察延迟与错误率再调整；
   不要无上限地并发打满（服务经由公网隧道，链路带宽有限）。
3. **复用连接**：使用支持 keep-alive 的 HTTP 客户端（如 `httpx`/`requests` 的 Session），
   避免每请求新建 TCP/TLS 连接。
4. **设置超时**：连接超时 5s，读取超时 30–60s（化学计算 + 跨公网链路，单请求可能 150–400ms 甚至更高）。
5. **重试要带退避**：仅对 `429/500/502/503/504` 和网络超时重试，采用指数退避 + 抖动
   （如 0.5s → 1s → 2s，最多 3–5 次）。对 `4xx`（`400/401/422`）**不要重试**，那是请求本身的问题。
6. **串行会很慢**：跨公网链路延迟较高，**请用并发**而非串行循环来提高吞吐。

---

## 4. 批量比较接口（推荐）

**请求** `POST /batch_mol_compare`

```json
{
  "pairs": [
    {"mol1": "CCO",     "mol2": "OCC",      "pair_id": "p1"},
    {"mol1": "aspirin", "mol2": "CC(=O)Oc1ccccc1C(=O)O", "pair_id": "p2"}
  ]
}
```
- `pair_id` 可选，用于在结果中对应每一对。

**响应**

```json
{
  "total": 2,
  "success": 2,
  "failed": 0,
  "results": [
    {"pair_id":"p1","mol1":"CCO","mol2":"OCC","exact_match":true,"tanimoto":1.0,"warning":null,"error":null},
    {"pair_id":"p2","mol1":"aspirin","mol2":"CC(=O)Oc1ccccc1C(=O)O","exact_match":true,"tanimoto":1.0,"warning":null,"error":null}
  ]
}
```
- 单对失败不会中断整体：失败项的 `error` 字段会给出原因，`failed` 计数加一。
- 单次请求的 `pairs` 数量建议 **≤ 200**；更多请分批提交。

---

## 5. 错误码对照

| HTTP | 含义 | 调用方处理 |
|------|------|-----------|
| 200 | 成功 | 仍需检查业务字段（如空 `smiles`、`results[].error`） |
| 401 | 认证失败 | 检查 `Authorization` 头与 key，**不重试** |
| 422 | 请求体字段/格式错误 | 修正请求，**不重试** |
| 500 | 服务内部错误（如解析失败） | 可少量退避重试；持续失败则检查输入 |
| 502 / 503 / 504 | 网关 / 不可用 / 超时 | 退避重试 |

---

## 6. 调用示例

### 6.1 Python（httpx，带并发限制 + 重试退避，推荐模板）

```python
import asyncio
import random
import httpx

BASE = "https://<服务域名或公网IP>/chemdraw/api"
API_KEY = "<API_KEY>"
HEADERS = {"Authorization": f"Bearer {API_KEY}", "Content-Type": "application/json"}

# 控制并发：按服务方建议设置（示例 30）
SEM = asyncio.Semaphore(30)
RETRIABLE = {429, 500, 502, 503, 504}


async def post(client: httpx.AsyncClient, path: str, payload: dict, max_retry: int = 4):
    async with SEM:
        for attempt in range(max_retry + 1):
            try:
                r = await client.post(f"{BASE}{path}", json=payload, headers=HEADERS)
                if r.status_code in RETRIABLE and attempt < max_retry:
                    raise httpx.HTTPStatusError("retriable", request=r.request, response=r)
                r.raise_for_status()
                return r.json()
            except (httpx.TransportError, httpx.HTTPStatusError) as e:
                # 4xx (except 429) should not be retried
                status = getattr(getattr(e, "response", None), "status_code", None)
                if status and status not in RETRIABLE:
                    raise
                if attempt == max_retry:
                    raise
                await asyncio.sleep((2 ** attempt) * 0.5 + random.random() * 0.3)


async def main():
    limits = httpx.Limits(max_connections=50, max_keepalive_connections=50)
    timeout = httpx.Timeout(connect=5.0, read=60.0, write=10.0, pool=5.0)
    async with httpx.AsyncClient(limits=limits, timeout=timeout) as client:
        # 单条：name -> smiles
        print(await post(client, "/name2smiles", {"name": "aspirin"}))

        # 批量比较（首选）：一次请求多对
        pairs = [{"mol1": "CCO", "mol2": "OCC", "pair_id": f"p{i}"} for i in range(50)]
        result = await post(client, "/batch_mol_compare", {"pairs": pairs})
        print(result["success"], "/", result["total"])


if __name__ == "__main__":
    asyncio.run(main())
```

### 6.2 curl

```bash
# 健康检查（无需认证）
curl -s https://<域名>/chemdraw/api/health

# 名称转 SMILES
curl -s -X POST https://<域名>/chemdraw/api/name2smiles \
  -H "Authorization: Bearer <API_KEY>" \
  -H "Content-Type: application/json" \
  -d '{"name": "aspirin"}'

# 批量比较
curl -s -X POST https://<域名>/chemdraw/api/batch_mol_compare \
  -H "Authorization: Bearer <API_KEY>" \
  -H "Content-Type: application/json" \
  -d '{"pairs":[{"mol1":"CCO","mol2":"OCC"}]}'
```

---

## 7. 给 Code Agent 的速记（TL;DR）

- 认证：`Authorization: Bearer <API_KEY>`；JSON 进 JSON 出。
- 多对比较 → **一律用 `batch_mol_compare`**，不要循环单调。
- 并发 20–50，连接 keep-alive，read 超时 60s。
- 重试只针对 5xx/超时，指数退避；4xx 直接报错。
- `name2smiles` 可能返回空字符串；`batch` 的每项可能带 `error`。
- 不要调用 `/smiles2rdkit`（暂不可用）。
