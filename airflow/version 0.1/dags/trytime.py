import time
import multiprocessing

for i in range(3):
    print(i)
    time.sleep(10)
    
print(multiprocessing.cpu_count())
