import pickle
import matplotlib.pyplot as plt
import seaborn as sns

# Load the UMAP embedding from the pickle file
def load_pickle(filename):
    return pickle.load(open(file=filename, mode='rb'))

outdir = '../results/crowdingDE'
umap_ordination = load_pickle(f"{outdir}/evo_umap.pkl")

# Extract embeddings (assuming umap_ordination is a tuple of (embedding, embedding))
embedding = umap_ordination[0]

# Plotting the UMAP results
plt.figure(figsize=(8, 6))
sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], s=50, alpha=0.7, edgecolor=None)
plt.title("UMAP Embedding")
plt.xlabel("UMAP Dimension 1")
plt.ylabel("UMAP Dimension 2")
plt.grid(True)
plt.savefig("../figures/test/umap_evo.png")
