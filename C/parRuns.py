import sys
import subprocess
import time

numOfRuns = int(sys.argv[1])
print numOfRuns

for i in range(numOfRuns):
    subprocess.call(['./Generator', 'initial'+str(i)+'.csv'])
    time.sleep(1)   #Let the random number generator seed change

for i in range(numOfRuns):
    subprocess.call(['mkdir', 'Results'+str(i)])
    subprocess.Popen(['gnome-terminal', '--', './LIF_2D_Classic', 'initial'+str(i)+'.csv', str(i)])
