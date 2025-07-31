# dependencies.py
from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from typing import Annotated
from config import settings

# 1. 创建一个 HTTPBearer 实例。它会自动查找 Bearer token。
security_scheme = HTTPBearer()

# 2. 按照您提供的模式重写验证函数
def verify_api_key(credentials: Annotated[HTTPAuthorizationCredentials, Depends(security_scheme)]):
    """
    验证 Bearer Token 格式的 API Key。
    这个函数会被用作一个依赖项来保护API路由。
    """
    # 如果 .env 文件中没有设置API_KEY，或者为默认值，则跳过验证（方便开发）
    if not settings.API_KEY or settings.API_KEY == "DEFAULT_API_KEY":
        return True

    # credentials.credentials 会包含 "Bearer " 后面的 token 部分
    if credentials.credentials != settings.API_KEY:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="无效的API Key",
            headers={"WWW-Authenticate": "Bearer"},
        )
    
    # 验证成功，返回 True。
    # 端点函数签名中的 `_: bool = Depends(...)` 会接收这个返回值，但通常会忽略它。
    return True
