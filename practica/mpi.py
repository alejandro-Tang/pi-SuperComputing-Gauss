import os
import sys
import matplotlib.pyplot as plt

def main():
    
    cores = int(input("Insert number of cores: "))
    clean_comand = 'make clean'
    compile_comand = 'make'
    os.system(clean_comand)
    os.system(compile_comand)
    test_time = []
    test_size = [8, 800, 1600, 2400, 3200]
    test_cores = []
    file_name = sys.argv[1]
    print("File name: ", file_name)

    for i in range(0, cores):
        out_time = os.popen('mpirun -np ' + str(2**i) + ' gaussian '+ file_name).read()
        test_cores.append(2**i)
        test_time.append(float(out_time[:-2]))
    print(test_time)
    plt.plot(test_cores,test_time)
    plt.title("Time taken for Gaussian Elimination.")
    plt.xlabel("Number of threads / AU", fontsize=10)
    plt.ylabel("Time / s", fontsize=10)
    for y in range(0, int(max(test_time))):
      plt.plot(range(1,2**(cores-1)+1), [y] * len(range(1,2**(cores-1)+1)), "--", lw=0.5, color="black", alpha=0.3)
    plt.xticks(range(1,2**(cores-1)+1), fontsize=10)  
    plt.yticks(fontsize=10)
    
    fig, ax = plt.subplots()
    ax.scatter(test_cores, test_time)

    for i, txt in enumerate(test_time):
      plt.annotate(txt, (test_cores[i], test_time[i]))
       
    plt.show()

    pause = input("Press any key to continue.")


main()

