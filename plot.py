import matplotlib.pyplot as plt
from Bio import Entrez
import pandas as pd
import time
import json
import os

# Function to fetch PubMed data
def fetch_pubmed_data(query, year_from, year_to, cache_filename='pubmed_cache.json'):
    if os.path.exists(cache_filename):
        with open(cache_filename, 'r') as file:
            cache = json.load(file)
    else:
        cache = {}

    results_per_year = {}
    for year in range(year_from, year_to + 1):
        cache_key = f"{query}_{year}"
        if cache_key in cache:
            results_per_year[year] = cache[cache_key]
        else:
            time.sleep(1)  # to avoid hitting API rate limits
            handle = Entrez.esearch(db="pubmed", term=query, mindate=year, maxdate=year, retmax=1000)
            record = Entrez.read(handle)
            handle.close()
            count = int(record['Count'])
            results_per_year[year] = count
            cache[cache_key] = count

    with open(cache_filename, 'w') as file:
        json.dump(cache, file)

    return results_per_year


# Fetch data for different Monte Carlo codes
year_from = 1989
year_to = 2024

geant4_data = fetch_pubmed_data("Geant4", year_from, year_to)
fluka_data = fetch_pubmed_data("FLUKA", year_from, year_to)
mcnp_data = fetch_pubmed_data("MCNP", year_from, year_to)
egsnrc_data = fetch_pubmed_data("EGSnrc", year_from, year_to)
penelope_data = fetch_pubmed_data("Penelope", year_from, year_to)

# Convert data into a Pandas DataFrame for easy plotting
years = list(range(year_from, year_to + 1))
data = pd.DataFrame({
    "Year": years,
    "Geant4": [geant4_data.get(year, 0) for year in years],
    "FLUKA": [fluka_data.get(year, 0) for year in years],
    "MCNP": [mcnp_data.get(year, 0) for year in years],
    "EGSnrc": [egsnrc_data.get(year, 0) for year in years],
    "Penelope": [penelope_data.get(year, 0) for year in years],
})

# Plotting the data with enhanced aesthetics
plt.figure(figsize=(10, 7))
plt.plot(data["Year"], data["Geant4"], '-o', label="Geant4", color='#1f77b4', linewidth=2, markersize=6, markerfacecolor='white', markeredgewidth=2)
plt.plot(data["Year"], data["FLUKA"], '-s', label="FLUKA", color='#ff7f0e', linewidth=2, markersize=6, markerfacecolor='white', markeredgewidth=2)
plt.plot(data["Year"], data["MCNP"], '-^', label="MCNP", color='#2ca02c', linewidth=2, markersize=6, markerfacecolor='white', markeredgewidth=2)
plt.plot(data["Year"], data["EGSnrc"], '-D', label="EGSnrc", color='#d62728', linewidth=2, markersize=6, markerfacecolor='white', markeredgewidth=2)
plt.plot(data["Year"], data["Penelope"], '-P', label="Penelope", color='#9467bd', linewidth=2, markersize=6, markerfacecolor='white', markeredgewidth=2)

# Enhancing the grid and ticks (less dense grid)
plt.minorticks_on()
#plt.grid(True, which='major', linestyle='--', linewidth=0.7, color='gray')
plt.grid(True, which='major', linestyle='--', linewidth=0.7, color='#cccccc') #light gray
plt.grid(True, which='minor', linestyle='', linewidth=0)  # Disable minor grid

# Customizing ticks to show fewer horizontal grid lines (two lines per tick label)
plt.yticks(range(0, 201, 20))

# Formatting the plot
plt.xlabel("Year", fontsize=16, fontweight='bold')
plt.ylabel("Number of Publications", fontsize=16, fontweight='bold')
plt.title("Publications related to Monte Carlo tools on PubMed", fontsize=18, fontweight='bold', pad=20)
plt.xlim(year_from, year_to)
plt.ylim(0, 200)
plt.xticks(range(year_from, year_to + 1, 5))

# Legend inside the plot
plt.legend(loc='upper left', fontsize=12, frameon=True, shadow=True, borderpad=1)

# Adjust layout to fit everything
plt.tight_layout()

# Save plot to file
plt.savefig('MC_PubMed.png')

# Show plot
#plt.show()
