# %%
import pandas as pd
import numpy as np
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# %%
# DIR = "/scratch/users/taerim/projects/YiLu_DNAzyme/v1_clean/results"
# sample = "aKG"
# round = 18
# seq_start = "TACGACTCACTATAGGG"
# seq_end = "ATCTGACGGTAACGCTATAGTGTCACCTAAAT"
# Deduplicate = "UseAll"
# Deduplicate = "Deduplicate"

DIR = sys.argv[1]
sample = sys.argv[2]
round = int(sys.argv[3])
seq_start = sys.argv[4]
seq_end = sys.argv[5]
Deduplicate = sys.argv[6]

# %%
DIR_Output = f"{DIR}/{sample}{round}/"
# print(DIR_Output)

# %%
seqs = pd.read_csv(
    # "/scratch/users/taerim/projects/YiLu_DNAzyme/v1/data/NGS_Nanopore_raw_data/aKG/aKG19_L001/tmp_strand12.txt",
    # "/scratch/users/taerim/projects/YiLu_DNAzyme/v1/results/aKG19/aKG19_PEAR.assembled.5to3.merged.final.txt",
    f"{DIR_Output}/{sample}{round}_PEAR.assembled.5to3.merged.final.txt",
    header=None,
)
seqs = seqs[0].tolist()
# print(len(seqs))
# seqs[:5]


# %%
def extract_subseq(seq, start, end):
    # if multiple start or end sequences are present, take the first one for start and the last one for end
    start_idx = seq.find(start)
    end_idx = seq.rfind(end)
    if start_idx != -1 and end_idx != -1 and start_idx < end_idx:
        return seq[start_idx : end_idx + len(end)]
    else:
        return None


# %%
seqs_clean = []
lens = []
for x in seqs:
    subseq = extract_subseq(x, seq_start, seq_end)
    if subseq is not None:
        seqs_clean.append(subseq)
        lens.append(len(subseq))

# %%
expected_length = pd.Series(lens).value_counts().index[0]
# print(f"Expected length: {expected_length}")

# %%
# print(len(seqs_clean))
seqs_clean_len = [s for s in seqs_clean if len(s) == expected_length]
# print(len(seqs_clean_len))

# %%
# Extract N40 sequences
seqs_clean_len_n40 = []
for s in seqs_clean_len:
    start_idx = s.find(seq_start)
    end_idx = s.rfind(seq_end)
    if start_idx != -1 and end_idx != -1 and start_idx < end_idx:
        n40_seq = s[start_idx + len(seq_start) : end_idx]
        seqs_clean_len_n40.append(n40_seq)

# %%
# Save N40 sequence counts.
if Deduplicate == "UseAll":
    pd.DataFrame(pd.Series(seqs_clean_len_n40).value_counts()).reset_index(
        drop=False
    ).rename(columns={"index": "N40Seq", 0: "counts"}).to_csv(
        f"{DIR_Output}/{sample}{round}_N40Seq_Counts_{Deduplicate}.tsv",
        index=False,
        sep="\t",
    )

# %%
# Save top N40 sequences as a bar plot.
if Deduplicate == "UseAll":
    Top_N = 20
    pd.Series(seqs_clean_len_n40).value_counts().head(Top_N).plot(
        kind="bar", figsize=(20, 10)
    )
    plt.title(f"Top {Top_N} N40 sequences for {sample} round {round}")
    plt.xlabel("N40 sequence")
    plt.ylabel("Counts")
    # rotate x-axis labels for better readability
    plt.xticks(rotation=15, ha="right", fontsize=7)

    # Add cululative frequency line (y range 0-100%)
    counts = pd.Series(seqs_clean_len_n40).value_counts().head(Top_N)
    cumulative_counts = counts.cumsum()
    cumulative_percentage = cumulative_counts / cumulative_counts.iloc[-1] * 100
    plt.twinx()
    plt.plot(
        cumulative_percentage, color="red", marker="o", label="Cumulative Frequency (%)"
    )
    plt.ylabel("Cumulative Frequency (%)")
    plt.ylim(0, 110)
    plt.legend(loc="upper right")

    plt.savefig(f"{DIR_Output}/{sample}{round}_N40Seq_Top{Top_N}.png")

# %%


# %%
df_counts = pd.concat(
    [
        pd.Series(seqs_clean_len_n40).value_counts().head(50),
        pd.Series(seqs_clean_len_n40).value_counts(normalize=True).head(50),
    ],
    axis=1,
    keys=["Counts", "Frequency"],
)
df_counts["Cumulative Frequency"] = df_counts["Frequency"].cumsum()
# df_counts.head(20)

# %%
df_seq = pd.DataFrame(
    {
        "sequence": pd.Series(seqs_clean_len_n40).value_counts().index,
        "counts": pd.Series(seqs_clean_len_n40).value_counts().values,
    }
)
# print(df_seq.shape)
# df_seq.head()

# %%
# Keep only sequences with counts above a certain threshold. It's an arbitry value, but it helps to reduce the number of sequences to analyze.
threshold = 10
df_seq_dedup = df_seq.loc[df_seq.counts >= threshold, :].copy()

if Deduplicate == "Deduplicate":
    print("Deduplicate N40Seq counts to 1 for each unique sequence.")
    df_seq_dedup["counts"] = 1

# print(df_seq_dedup.shape)
# df_seq_dedup.head()

# %%
print(df_seq_dedup.shape)
df_seq_dedup.head()

# %% [markdown]
# # Make 2D matrix (motif x position)

# %%
k = 3
bases = "ACGT"

motifs = ["".join(p) for p in itertools.product(bases, repeat=k)]
# print(f"Total {k}-mer motifs: {len(motifs)}")
# print(motifs[:5])

# %%


# %%
# Make template DataFrame for motif counting
motif_position = []
for m in motifs:
    for p in range(len(df_seq_dedup.sequence[0]) - k + 1):
        motif_position.append(f"{m}_{p}")

df_template = pd.DataFrame({"motif_position": motif_position})
df_template["motif"] = df_template.motif_position.apply(lambda x: x.split("_")[0])
df_template["position"] = df_template.motif_position.apply(
    lambda x: int(x.split("_")[1])
)
# df_template["value"] = 0

# df_template.head()

# %%
# count motif_position
motif_position_count = []

for seq_tmp, count in zip(df_seq_dedup.sequence, df_seq_dedup.counts):
    motif_position_count_tmp = []
    for i in range(len(seq_tmp) - k + 1):
        kmer = seq_tmp[i : i + k]
        motif_position_count_tmp.append(f"{kmer}_{i}")
        motif_position_count_tmp_df = pd.DataFrame(
            {
                "motif_position": motif_position_count_tmp,
                "count": count,
            }
        )
    motif_position_count.append(motif_position_count_tmp_df)

motif_position_count_concat = pd.concat(motif_position_count, ignore_index=True)
# print(motif_position_count_concat.shape)
# motif_position_count_concat.head()

# %%
motif_position_count_concat = motif_position_count_concat.groupby(
    "motif_position"
).sum()
motif_position_count_concat.sort_values("count", ascending=False, inplace=True)
# motif_position_count_concat.head()

# %%


# %%
df_final = df_template.merge(
    motif_position_count_concat,
    left_on="motif_position",
    right_on="motif_position",
    how="left",
)
df_final["count"] = df_final["count"].fillna(0).astype(int)
# df_final.head()

# %%
# write motif/count table.
df_final.to_csv(
    f"{DIR_Output}/{sample}{round}_N40Seq_MotifPositionCount_{Deduplicate}.tsv",
    index=False,
    sep="\t",
)

# %%
df_final_2d = df_final.pivot(index="motif", columns="position", values="count")
df_final_2d.head()

# %%


# %%
# 2D heatmap
plt.figure(figsize=(12, 12))
sns.heatmap(
    df_final_2d / df_final_2d.sum(axis=0),  # normalize by column sum
    cmap="Reds",
    vmin=0,
    vmax=1,
)
plt.title(f"2D Heatmap for {sample} round {round}")
plt.xlabel("Position")
plt.ylabel("Motif")
plt.savefig(
    f"{DIR_Output}/{sample}{round}_N40Seq_MotifPosition_2DHeatMap_{Deduplicate}.png"
)
