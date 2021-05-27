## 本地查询基因的标准名称

### 程序设计思想

#### · 数据库处理

    - 利用官网下载的数据库文件，构建本地的json数据库（好像有个pickle更节省存储）
    - json，节省每次查询的字典生成时间

#### · 数据库结构
    - 如果status为symbol withdrawn，则舍弃改行；因为该描述已经合并到其他结果中
    - 如果一个基因，存在previous_symbol，该previous同样需要作为key

#### · 查询输入
    1. 单个基因
    2. 多个基因（竖线或者逗号分隔）
    3. 输入基因列表

#### · 结果输出
    1. 直接本地展示输出
    2. 导出结果至excel表格中
    3. 使用pd.DataFrame展示输出结果
        
    
| query_gene      | HGNC | HGNC_previous | status |
| ----------- | ----------- | ------- | --------- |
| EGFR      | EGFR       | ERBB | Approved |
| CORTBP1   | SHANK2        | CORTBP1 | Approved |