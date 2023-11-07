import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from Nal5 import corelation_FFT


with open("Naloga5\\significant-volcanic-eruption-database.csv", mode="r") as f:
    lines = f.readlines()
    eruptions = dict()
    for line in lines[1:]:
        year = int(line.split(";")[0])

        try:
            eruptions[year] = eruptions[year] + 1
        except KeyError:
            eruptions[year] = 1

    
years = np.arange(min(eruptions.keys()), max(eruptions.keys()))
years = np.arange(1800, max(eruptions.keys()))

erupts = []
for year in years:
    try:
        erupts.append(eruptions[year])
    except KeyError:
        erupts.append(0)

plt.plot(years, erupts)
plt.xlabel("Leto")
plt.ylabel("Število izbruhov v letu")
plt.show()

ns, corr =  corelation_FFT(erupts,erupts)

plt.plot(ns, corr)
plt.xlabel("Zamik")
plt.show()




import pickle

"""
### ROČNA KLASIFIKACIJA DRŽAV IN REGIJ
regions = dict()
with open("Naloga5\\significant-volcanic-eruption-database.csv", mode="r") as f:
    lines = f.readlines()
    eruptions = dict()
    for line in lines[1:]:
        country = line.split(";")[7]

        if country not in regions.keys():
            regions[country] = input(country + "\n")

print(regions)

with open('Naloga5\\regions of countries.pkl', 'wb') as fp:
    pickle.dump(regions, fp)
    print('dictionary saved successfully to file')
"""

with open('Naloga5\\regions of countries.pkl', 'rb') as fp:
    regions = pickle.load(fp)

with open("Naloga5\\significant-volcanic-eruption-database.csv", mode="r") as f:
    lines = f.readlines()
    eruptionsNA = dict()
    eruptionsSA = dict()
    eruptionsEU = dict()
    eruptionsAS = dict()
    eruptionsAF = dict()
    eruptionsAU = dict()
    for line in lines[1:]:
        year = int(line.split(";")[0])
        country = line.split(";")[7]

        if regions[country] == "NA":
            try:
                eruptionsNA[year] = eruptionsNA[year] + 1
            except KeyError:
                eruptionsNA[year] = 1

        if regions[country] == "EU":
            try:
                eruptionsEU[year] = eruptionsEU[year] + 1
            except KeyError:
                eruptionsEU[year] = 1

        if regions[country] == "AS":
            try:
                eruptionsAS[year] = eruptionsAS[year] + 1
            except KeyError:
                eruptionsAS[year] = 1

        if regions[country] == "AF":
            try:
                eruptionsAF[year] = eruptionsAF[year] + 1
            except KeyError:
                eruptionsAF[year] = 1

        if regions[country] == "SA":
            try:
                eruptionsSA[year] = eruptionsSA[year] + 1
            except KeyError:
                eruptionsSA[year] = 1

    eruptions = [eruptionsSA, eruptionsAF, eruptionsAS, eruptionsEU, eruptionsNA]
    titles = ["Južna Amerika", "Afrika", "Azija", "Evropa", "Severna Amerika"]

for i in range(len(eruptions)):
    years = np.arange(1900, 2024)

    eruptsI = []
    for year in years:
        try:
            eruptsI.append(eruptions[i][year])
        except KeyError:
            eruptsI.append(0)

    for j in range(len(eruptions)):    
        eruptsJ = []
        for year in years:
            try:
                eruptsJ.append(eruptions[j][year])
            except KeyError:
                eruptsJ.append(0)

        #plt.plot(years, eruptsI)
        #plt.plot(years, eruptsJ)
        #plt.title(f"{titles[i]} {titles[j]}")
        #plt.show()

        ns, corr =  corelation_FFT(eruptsI,eruptsJ)
        corr = corr / (np.sum(eruptsI) + np.sum(eruptsJ))

        plt.plot(ns, corr, label = f"Korelacija {titles[i]} {titles[j]}")
        plt.xlabel("Zamik")
    plt.legend()
    plt.show()

    