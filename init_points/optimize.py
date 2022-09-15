import numpy as np
num_points_in_equation = 7


with open("/home/behnam/phd/Research/ICP/init points/points_backup") as f:
	lines = f.readlines()
ref = []
query = []
for i in range(0,len(lines),6):
	temp_ref = []
	temp_query = []
	for j in range(3):
		temp_ref.append(lines[i+j].split("\n")[0])
		temp_query.append(lines[i+j+3].split("\n")[0])
	ref.append(temp_ref)
	query.append(temp_query)
ref = np.array(ref).astype(np.float64)
query = np.array(query).astype(np.float64)

#distances = np.sqrt(((ref - query)^2))
ref = ref[0:num_points_in_equation]
query = query[0:num_points_in_equation]


distances = np.sum(np.sqrt((pow((ref - query),2))),1)





A = ref - query
#b = np.zeros((num_points_in_equation, 1))






x = np.linalg.lstsq(A,distances)
#x = np.linalg.solve(A,b)


