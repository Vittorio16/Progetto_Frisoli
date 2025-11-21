g = 9.81
# Questi parametri possono essere forniti a run time, importante mantenerli in formato espandibile
NUM_PENDOLI = 2
MASSE = [1,1]
LUNGHEZZE = [1,1]
A = 5
w = 2

# Stato iniziale - z = [theta_1, ......, theta_n, d(theta_1)/dt, ......., d(theta_n)/dt] 
# Angoli in gradi
state_vector = [5,10,2,3]

