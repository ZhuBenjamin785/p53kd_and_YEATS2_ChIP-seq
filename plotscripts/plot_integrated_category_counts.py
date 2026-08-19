"""Plot corrected integrated category counts."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plot_common import BASE,save
def main():
    d=pd.read_csv(BASE/"integration/tables/significant_category_counts.csv"); fig,ax=plt.subplots(figsize=(10,6)); sns.barplot(data=d,x="category",y="count",hue="scope",ax=ax); ax.set_title("Integrated significant category counts"); ax.tick_params(axis="x",rotation=25); save(fig,"integrated_category_counts")
if __name__=="__main__": main()
