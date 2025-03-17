import numpy as np

ages0_3 = [30,40,55,80,70,30,60,90,60,35,70,60,44,89,52,40,28,60,72,68]
times0_3 = [10, 17, 5, 7, 30, 12, 5, 22, 25, 9, 25, 15, 11, 35, 8, 15, 10, 27, 12, 30]
times0_4 = [10,20,8,12,32,16,11,21,36,14,28,24,18,28,13,24,16,36,45,56]
ages0_4 = [30,40,55,74,70,30,60,72,60,35,75,60,44,65,52,40,28,60,72,68]
times0_5 = [10,25,10,15,40,20,10,26,45,18,35,30,23,36,16,31,20,45,25,56]
ages0_5 = [30,40,55,74,70,30,60,72,60,35,75,60,44,65,52,40,28,60,72,68]
times0_6 = [20,35,10,15,60,25,10,44,50,18,50,30,23,70,16,31,20,55,25,61]
ages0_6 = [30,40,55,80,70,30,60,90,60,35,70,60,44,89,52,40,28,60,72,68]
times0_7 = [23,36,12,17,70,29,12,51,58,21,58,35,27,81,19,36,23,46,29,64]
ages0_7 = [30,40,55,80,70,30,60,90,60,35,70,60,44,89,52,40,28,60,72,68]

print(np.std((np.array(ages0_3)- np.array(times0_3))))
print(np.std(ages0_3))
print(np.std(times0_3))

print(np.std((np.array(ages0_4)- np.array(times0_4))))
print(np.std(ages0_4))
print(np.std(times0_4))
      
print(np.std((np.array(ages0_5)- np.array(times0_5))))
print(np.std(ages0_5))
print(np.std(times0_5))

print(np.std((np.array(ages0_6)- np.array(times0_6))))
print(np.std(ages0_6))
print(np.std(times0_6))
      
print(np.std((np.array(ages0_7)- np.array(times0_7))))
print(np.std(ages0_7))
print(np.std(times0_7))

print('std_event_times')
print(np.std(times0_3))
print(np.std(times0_4))
print(np.std(times0_5))
print(np.std(times0_6))
print(np.std(times0_7))

print('mean_evnt_times')

print(np.mean(times0_3))
print(np.mean(times0_4))
print(np.mean(times0_5))
print(np.mean(times0_6))
print(np.mean(times0_7))
