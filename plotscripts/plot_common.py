from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)
sns.set_theme(style="whitegrid", context="talk")

def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)

def note(name, message):
    print(f"[{name}] {message}")
