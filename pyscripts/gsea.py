
# %%
import pandas as pd
import mygene
import pydeseq2.preprocessing

# %%
# metadata = pd.read_csv("../HCC/SraRunTable(9).csv", set="\t").filter(items=["Run", "tissue"])
data = pd.read_csv("../deseq2_norm_expression.gct")
data = data.set_index(data.iloc[:, 0])#.drop("Unnamed: 0", axis=1)


def expression_matrix_prep(data, datatype):
    mg = mygene.MyGeneInfo()
    annotated_data = data.copy()

    # Map full gene ID → base ID (stripped of .xx)
    # full_to_base = {gene for gene in data.index} and then u split yah

    # results = mg.querymany(full_to_base.values(), fields="symbol", species="human")
    results = mg.querymany(
        annotated_data.index.to_list(), fields="symbol", species="human"
    )

    # mapping the gene id to symbol
    id_to_symbol = {}
    not_annotated = []
    for r in results:
        if r.get("notfound"):
            not_annotated.append(r["query"])
            id_to_symbol[r["query"]] = "NONE"
        elif "symbol" in r:
            id_to_symbol[r["query"]] = r["symbol"]

    # Assign annotations efficiently
    annotated_data.loc[data.index, "gene_symbol"] = [
        id_to_symbol.get(gene, None) for gene in data.index
    ]

    # Append not-annotated gene IDs to file
    with open(("../%s/non_annotated.txt" % datatype), "a") as f:
        for gene in data.index:
            if gene in not_annotated:
                f.write(gene + "\n")

    column_to_move = "gene_symbol"
    column_value = annotated_data.pop(column_to_move)
    annotated_data.insert(0, column_to_move, column_value)
    annotated_data.reset_index(inplace=True)
    annotated_data.set_index("gene_symbol", inplace=True)
    return annotated_data


# %%
df = expression_matrix_prep(data, "HCC")
df.rename(columns={"Unnamed: 0": "gene_id"}, inplace=True)
# %% dropping ones without gene annotation
df = df[(df.index.notna()) & (df.index != "None") & (df.index != "NONE")]
df = df[(df != 0).any(axis=1)]
# %%
df_norm, size_factors = pydeseq2.preprocessing.deseq2_norm(
    df.loc[:, "SRR3140234":"SRR3140340"].T
)
df.loc[:, "SRR3140234":"SRR3140340"] = df_norm.T
# %%
# col_b = df.pop("SRR3140234")
# df.insert(len(df.columns), "SRR3140234", col_b)
df.to_csv("../HCC/expression_matrix.gct", sep="\t")

# %%
### ===== GSEAPY ===== ###
