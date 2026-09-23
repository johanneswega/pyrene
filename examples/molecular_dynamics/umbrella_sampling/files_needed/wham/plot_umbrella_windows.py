from pyrene.standard.packages import *
from pyrene.standard.misc import rainbow

# ============================================================
# SETTINGS
# ============================================================

HISTOGRAM_FILE = "histogram.xvg"

# Window centers, in nm.
# Update this list to match the windows included in your WHAM run.
WINDOW_CENTERS = np.arange(0.35, 2.05, 0.05)


# ============================================================
# READ XVG DATA
# ============================================================

data = np.loadtxt(
    HISTOGRAM_FILE,
    comments=("@", "#")
)

r = data[:, 0]
histograms = data[:, 1:]

print(f"Number of coordinate bins: {len(r)}")
print(f"Number of histogram columns: {histograms.shape[1]}")
print(f"Number of supplied window centers: {len(WINDOW_CENTERS)}")


# ============================================================
# CHECK CONSISTENCY
# ============================================================

if histograms.shape[1] != len(WINDOW_CENTERS):
    raise ValueError(
        "The number of histogram columns does not match "
        "the number of window centers."
    )


# ============================================================
# PLOT ALL HISTOGRAMS
# ============================================================

fig, ax = plt.subplots(figsize=(10, 6))

colors = rainbow(WINDOW_CENTERS)

for i, center in enumerate(WINDOW_CENTERS):
    ax.plot(
        r,
        histograms[:, i],
        linewidth=1.2,
        label=f"{center:.2f} nm",
        color=colors[i]
    )

ax.set_xlabel("BPE–DPA COM distance / nm")
ax.set_ylabel("$N$")
ax.set_title("Umbrella sampling histograms")

ax.legend(
    title="Window",
    fontsize=8,
    ncol=2
)

#ax.grid(True, alpha=0.3)
fig.tight_layout()

plt.savefig("histograms_all.svg", transparent=True)
plt.show()