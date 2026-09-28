# MicrobiomeAnalyst Prep

## Interface gráfica

A aplicação disponibiliza um preparador genérico de dados para a Fase 2. Não existem
grupos, organismos ou nomes de amostras definidos no código da interface.

O utilizador carrega:
- uma tabela de abundâncias (`.csv`, `.tsv`, `.xls`, `.xlsx` ou `.xlsm`);
- uma tabela de taxonomia nos mesmos formatos;
- um ficheiro FASTA.

Depois pode associar as amostras a grupos de duas formas:
1. **Regras por prefixo** — uma regra é aplicada aos nomes que começam pelo prefixo indicado;
2. **Atribuição direta por amostra** — cada amostra pode receber um grupo individualmente.

Em ambos os modos existe um botão **"Sugerir a partir dos nomes das amostras"**
que propõe grupos automaticamente, olhando só para os nomes de amostra que
vêm do ficheiro carregado (não há nenhum grupo pré-definido no código). A
sugestão é sempre editável antes de avançar.

O resultado é composto por:
- `otu_table.txt`
- `otu_table_tss.txt`
- `taxonomy.txt`
- `metadata.txt`

e pode ser descarregado num único ZIP.

## Instalação

```bash
pip install -r requirements.txt
```

## Executar a interface

```bash
streamlit run app.py
```

A aplicação será aberta no navegador.

## Execução por linha de comandos

A funcionalidade original continua disponível através de:

```bash
python scripts/prep.py \
  --input data/feature-table.tsv \
  --fasta data/feature.fasta \
  --taxonomy data/taxonomy.tsv \
  --output outputs
```

O `config/mapping.json` pode continuar a ser usado para o modo de linha de comandos.
A interface gráfica não depende desse ficheiro.
