# %%
import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import umap

sns.set(style='white', context='poster', rc={'figure.figsize':(14,10)})

data = pd.read_table("/media/B_16TB/jason/CPE_analysis/test/VF_matrix.tsv", header=0)
data.head()

reducer = umap.UMAP()

# %%
scaled_data = StandardScaler().fit_transform(data.iloc[:, 2:])

embedding = reducer.fit_transform(scaled_data, n_neighbors=15, min_dist=0.8, metric='jaccard')
embedding.shape

plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=pd.Categorical(data["Source"]).codes,
    cmap='tab20'  # Use a larger colormap for more unique colors
)

# Add legend for color groups, placed outside the plot
handles = []
categories = pd.Categorical(data["Source"]).categories
for i, category in enumerate(categories):
    handles.append(plt.Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=plt.cm.tab20(i % 20), label=category, markersize=10))
plt.legend(handles=handles, title="Source", loc='center left', bbox_to_anchor=(1, 0.5))

plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP projection of the test dataset', fontsize=24)
plt.tight_layout()  # Adjust layout to make room for legend
plt.show()

