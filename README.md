# ChemDraw Server

> English | [中文](#中文说明)

## Overview

ChemDraw Server is a unified chemical informatics API service based on FastAPI and FastMCP. It provides endpoints for chemical name/SMILES conversion and molecule structure processing, suitable for integration into chemical drawing tools, automation pipelines, or as a backend for chemical informatics applications.

## Features

- Convert chemical names to SMILES strings
- Convert SMILES strings to chemical names
- Convert SMILES to RDKit molecule objects
- MCP (Multi-Chain Protocol) compatible API
- API Key authentication
- Prometheus monitoring integration

## Requirements

- Python == 3.10

## Installation

```bash
git clone https://github.com/tom832/chemdraw-server.git
cd chemdraw-server
uv sync
```

## Usage

### Start the server

```bash
uv run main_server.py
```

The API will be available at: `http://localhost:1145/chemdraw/api/`

### API Endpoints

- `POST /chemdraw/api/name_to_smiles`  
  Convert chemical name to SMILES

- `POST /chemdraw/api/smiles_to_name`  
  Convert SMILES to chemical name

- `POST /chemdraw/api/smiles_to_rdkit`  
  Convert SMILES to RDKit molecule

- `GET /chemdraw/api/health`  
  Health check

### API Key

Set your API key in a `.env` file or via environment variable:

```
API_KEY=your_api_key_here
```

## Configuration

See `config.py` for all configurable options.

## Dependencies

- fastapi
- fastmcp
- loguru
- mcp[cli]
- prometheus-fastapi-instrumentator
- rdkit
- uvicorn-loguru-integration

## License

MIT

---

## 中文说明

> [English](#overview) | 中文

### 项目简介

ChemDraw Server 是一个基于 FastAPI 和 FastMCP 的统一化学信息学 API 服务，提供化学名称/SMILES 互转、分子结构处理等接口，适用于集成到化学绘图工具、自动化流程或作为化学信息学应用的后端。

### 功能特性

- 化学名称转 SMILES
- SMILES 转化学名称
- SMILES 转 RDKit 分子对象
- 兼容 MCP 协议的 API
- API Key 认证
- 集成 Prometheus 监控

### 环境要求

- Python == 3.10

### 安装方法

```bash
git clone https://github.com/tom832/chemdraw-server.git
cd chemdraw-server
uv sync
```

### 启动服务

```bash
uv run main_server.py
```

API 默认地址为：`http://localhost:1145/chemdraw/api/`

### 主要接口

- `POST /chemdraw/api/name_to_smiles`  
  化学名称转 SMILES

- `POST /chemdraw/api/smiles_to_name`  
  SMILES 转化学名称

- `POST /chemdraw/api/smiles_to_rdkit`  
  SMILES 转 RDKit 分子对象

- `GET /chemdraw/api/health`  
  健康检查

### API Key 配置

在 `.env` 文件或环境变量中设置：

```
API_KEY=your_api_key_here
```

### 配置项

详见 `config.py` 文件。

### 依赖列表

- fastapi
- fastmcp
- loguru
- mcp[cli]
- prometheus-fastapi-instrumentator
- rdkit
- uvicorn-loguru-integration

### 许可证

MIT

---

如需进一步完善或定制内容，请告知！

---

# ⚠️ Copyright & Disclaimer | 版权与免责声明

**English:**

ChemDraw, ChemScript and related software are copyrighted by their respective companies. Please ensure you have purchased and are legally using genuine software. This project is for learning and academic reference only, and must not be used for any commercial purpose.

**中文：**

ChemDraw、ChemScript 等相关软件的版权归其所属公司所有。请确保您已购买并合法使用正版软件。本项目仅供学习与学术参考，不得用于任何商业用途。
