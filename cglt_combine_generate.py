import csv
import random
import numpy as np
from sklearn.cluster import KMeans

def read_combinations(filename='valid_combinations.csv'):
    combinations = []
    with open(filename, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            combinations.append((
                int(row['Sectionalradius']),
                int(row['angle']),
                int(row['radius']),
                int(row['Lumbus'])
            ))
    return combinations

def select_representative_pairs(combinations, num_pairs=10000):
    
    data = np.array(combinations)
    total_combinations = len(data)

    
    num_clusters = min(total_combinations, 500)  
    kmeans = KMeans(n_clusters=num_clusters, random_state=42).fit(data)
    cluster_centers = kmeans.cluster_centers_

    
    pairs = []
    while len(pairs) < num_pairs:
        
        cluster_indices = random.sample(range(num_clusters), min(num_clusters, 50))
        for cluster_idx in cluster_indices:
            combination1 = data[random.choice(np.where(kmeans.labels_ == cluster_idx)[0])]
            combination2 = data[random.randint(0, total_combinations - 1)]  
            jobnum = len(pairs) + 1
            pairs.append((*combination1, *combination2, jobnum))

            if len(pairs) >= num_pairs:
                break

    return pairs

def save_pairs_to_csv(pairs, filename='selected_pairs.csv'):
    with open(filename, mode='w', newline='') as file:  
        writer = csv.writer(file)
        
        writer.writerow([
            'Sectionalradius1', 'angle1', 'radius1', 'Lumbus1',
            'Sectionalradius2', 'angle2', 'radius2', 'Lumbus2',
            'jobnum'
        ])
        
        writer.writerows(pairs)


if __name__ == "__main__":
 
    combinations = read_combinations()

  
    num_pairs = 10000 
    pairs = select_representative_pairs(combinations, num_pairs)

    save_pairs_to_csv(pairs)
    print("已保存 {} 对组合到 selected_pairs.csv".format(len(pairs)))
