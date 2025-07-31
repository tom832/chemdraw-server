# config.py
from pydantic_settings import BaseSettings, SettingsConfigDict

class Settings(BaseSettings):
    # model_config 用于指定从 .env 文件加载环境变量
    model_config = SettingsConfigDict(env_file=".env", env_file_encoding='utf-8')
    
    # 定义环境变量
    API_KEY: str = "DEFAULT_API_KEY"
    PROJECT_NAME: str = "Chemdraw Tools API"
    PROJECT_VERSION: str = "1.0.0"

# 创建一个全局可用的配置实例
settings = Settings()