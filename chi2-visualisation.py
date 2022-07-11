import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

# Paramètres :
n = 8492
n2 = int(n / 2)
nb_tot = int((n * (n - 1)) / 2)

# Lecture des données (statistiques du chi-2)
# output_tot = open("Data/chi2-visualisation/output_tot_norejection_100.txt", 'r')
# output_matched = open("Data/chi2-visualisation/output_matched_norejection_100.txt", 'r')

output_tot = open("Data/chi2-visualisation/output_tot_reference_100.txt", 'r')
output_matched = open("Data/chi2-visualisation/output_matched_reference_100.txt", 'r')

# Premier graphe
t_tot = []
for ligne_tot in output_tot:
    temp_tot = float(ligne_tot)
    t_tot.append(temp_tot)

t_matched = []
for ligne_matched in output_matched:
    temp_matched = float(ligne_matched)
    t_matched.append(temp_matched)

X_tot = [i/nb_tot for i in range(nb_tot)]
X_matched = [i/n2 for i in range(n2)]

t_tot.sort()
t_matched.sort()

X_theo = np.arange(0, 100, 0.001)

fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(X_theo, chi2.cdf(X_theo, df=8), color = "green")
ax1.plot(t_tot, X_tot, color="red")
ax1.plot(t_matched, X_matched, color='blue')

# Deuxième graphe
t_matched_inv = dict()
for i in range(len(t_matched)):
    if (t_matched[i]  not in t_matched_inv):
        t_matched_inv[t_matched[i]] = [i]
    else:
        t_matched_inv[t_matched[i]].append(i)
    # print(i)

L = []
t_matched_set = set(t_matched)

count = 0
index = dict()
for t in t_matched:
    index[t] = 0

for t in t_tot:
    if t in t_matched_set and index[t] < len(t_matched_inv[t]):
        L.append((t,  t_matched_inv[t][index[t]]))
        index[t] += 1
    else:
        L.append((t, -1)) 

L2 = [0] * len(t_matched)
for i in range(len(L)):
    x, y = L[i]
    if y >= 0:
        L2[y] = i

X = [i for i in range(len(L2))]
ax2.scatter(X, L2)
ax2.plot([0.0, float(n2)], [0.0, float(nb_tot)], 'r-', lw=2)

plt.show()