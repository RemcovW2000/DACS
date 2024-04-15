dictt = {'1': 100, '2': 1292, '3': 88}  
  
# Getting the key with maximum value  
key_max = max(zip(dictt.values(), dictt.keys()))[1]  
value_max = max(zip(dictt.values(), dictt.keys()))[0] 

#value2 = max(dict.items(), key =lambda x:x[1])


print(key_max, value_max)
