import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import random
#fitting function
from scipy.optimize import curve_fit
from matplotlib.pyplot import imshow, show, colorbar







#gaussian function for fitting 
def gaussian(t,amp,mean,std):
    return amp*np.exp(-1/2*((t-mean)/std)**2)





time_res=100e-12  #50ps
total_time=2e-3 # unit: s
mesh_num=int(total_time/time_res)
print(mesh_num)
detector_tmng_res=350e-12
detector_efficiency=0.68

detector_dead_time=30e-9

single_photon_rate=10e+6*total_time #2.4e+6*total_time #unit: /s*s= number of photon events 

herald_efficiency=0.6 # this herald efficiency excludes Alice's APD detection efficiency 

dark_count_rate=500*total_time # unit: /s * s = number of photon events 

bg_count_rate=900e+3*total_time #900e+3 # unit: /s * s = number of photon events this is from atmosphere 

bg_counts_rate_from_unpair_photons=int(single_photon_rate*(1-herald_efficiency))



apd_dead_time=22e-9#84e-9 #unit of second

apd_dead_mesh=int(apd_dead_time/time_res)

time_tag_jitter=8e-12 #42e-12 #unit of second 

APD_jitter=350e-12 #unit of second 

system_jitters=np.sqrt(time_tag_jitter**2+APD_jitter**2)

total_system_jitters=np.sqrt(2)*system_jitters # this is the total system jitter from Ailice and Bob 



fileter_loss=0.9 # Filter loss

detection_efficiency=0.68

collecion_efficiency=0.5 #telescope, this is for Bob

atms_atten=-18 # unit dB and negative numbeer

channel_loss= 10**(atms_atten/10) #=0.0364 approxs -14 dB loss for downlink with 318 mm telescope  # this is for bob from the atmoshperic attenuation



total_efficiency=fileter_loss*detection_efficiency*herald_efficiency

bob_efficiency=channel_loss*collecion_efficiency





def random_num_gen(num):

  random_array=[]

  for i in range(num):

    random_array.append(random.randint(0,mesh_num))

  return random_array









photon_counts=np.zeros(mesh_num,dtype=np.uint8)



#random_dark_counts=random_num_gen(int(dark_count_rate)) #for alice and bob







random_single_photon_counts=np.random.randint(mesh_num,size=int(single_photon_rate)) #random_num_gen(int(single_photon_rate))





def loss_of_photon(array,efficiency):  #choose a list of photons with the efficiency 

    survived_photon=int(len(array)*efficiency)

    survived_array=random.sample(array,survived_photon)

    return survived_array



random_single_photon_counts_loss=loss_of_photon(list(random_single_photon_counts), total_efficiency)

random_single_photon_counts_loss_bob=np.copy(random_single_photon_counts_loss)

random_single_photon_counts_loss_bob=list(random_single_photon_counts_loss_bob)

random_single_photon_counts_loss_Bob=loss_of_photon(random_single_photon_counts_loss_bob, bob_efficiency)







alice_photon_counts= np.copy(photon_counts)

for i in range(0,len(random_single_photon_counts_loss)):

  alice_photon_counts[random_single_photon_counts_loss]=1

  

bob_photon_counts= np.copy(photon_counts)



for i in range(0,len(random_single_photon_counts_loss)):

  bob_photon_counts[random_single_photon_counts_loss_Bob]=1







# add back ground counts to bob's array 



random_bg_counts=np.random.randint(mesh_num,size=int(bg_count_rate+dark_count_rate+bg_counts_rate_from_unpair_photons))  #random_num_gen(int(bg_count_rate+dark_count_rate+bg_counts_rate_from_unpair_photons))  #for bob



for i in range(0,len(random_bg_counts)):

  bob_photon_counts[random_bg_counts]=1



#shift bob's array with respect to Alice's

  

time_shift=900e-12 #unit of second (s)  #shift bobs array by time_shift

shift_mesh=int(time_shift/time_res)





bob_photon_counts_roll=np.roll(bob_photon_counts,shift_mesh)

bob_photon_counts_roll[:shift_mesh]=0



alice_photon_counts_apd_dead=np.copy(alice_photon_counts)



#consider a dead time of APD 

for i in range(0,len(alice_photon_counts_apd_dead)):

    if alice_photon_counts_apd_dead[i]==1:

        alice_photon_counts_apd_dead[i+1:i+1+apd_dead_mesh]=0

    if bob_photon_counts_roll[i]==1:

        bob_photon_counts_roll[i+1:i+1+apd_dead_mesh]=0









bob_signal_index=[]

for i in range(0,len(bob_photon_counts_roll)):

    if int(bob_photon_counts_roll[i])==1:

        bob_signal_index.append(i)



def system_jitter(photon_timing):

    system_time_jitter=photon_timing+np.random.normal(0,total_system_jitters) #random time jitter due to the APD and time tagger

    return system_time_jitter

    





b_index=[]

for i in range(0,len(bob_photon_counts_roll)):

    if bob_photon_counts_roll[i]==1:

        b_index.append(i)





a_photon=[]





n=6 #how many adjecent time stamp we want to record (n-indexes preceding and trailing each Bob time)





for i in b_index:

    k=1

    total=0

    while total <n and i+k<len(alice_photon_counts):

        if int(alice_photon_counts[k+i])==1:

            a_photon.append(k)

            total+=1

        k+=1

        

    k=1

    total=0

    while total <n and i-k>0:

        if int(alice_photon_counts[i-k])==1:

            a_photon.append(-k)

            total+=1

        k+=1



a_photon_real_time=np.array(
    )*time_res



a_photon_with_time_jitter=np.zeros(len(a_photon_real_time))

for i in range(0,len(a_photon)):

    a_photon_with_time_jitter[i]=system_jitter(a_photon_real_time[i])











bin_resolution=10 #meaning that how much fraction 1 ns will be bin_size  

bins=np.array([i for i in range(-10*bin_resolution,10*bin_resolution)])*1/bin_resolution*1e-9

hist_data=plt.hist(a_photon_with_time_jitter,bins=bins,color='black',alpha=0)



plt.xticks(np.array([-10,-5,-2.5,-0.9,0,2.5,5,10])*1e-9,np.array([-10,-5,-2.5,-0.9,0,2.5,5,10]))

x_data=hist_data[1][0:-1]+time_res/2

y_data=hist_data[0]

print('number of hist_data={}'.format(len(x_data)))





# normalise the nosie offset 

# set the windows that are defined as the noise (from +/-2 ps are defined as nosie photon counts )

total_bg_counts=np.sum(y_data[int(len(y_data)/2):])*2

avg_bg_counts=total_bg_counts/len(y_data)





y_data_bg_corrected=y_data-avg_bg_counts





param_init=[np.max(y_data_bg_corrected),-1*time_shift,-1*total_system_jitters]



popt,pocv=curve_fit(gaussian,x_data,y_data_bg_corrected,p0=param_init)



x_=np.linspace(-10e-9,10e-9,1000)



plt.plot([popt[1],popt[1]],[0,popt[0]],color='red',linestyle='dotted')

plt.legend()

print('fitted_amplitude= {}counts'.format(popt[0]))

print('time difference = {}s'.format(popt[1]))

print('time resolution = {}s'.format(popt[2]))

print('n(adjacent photons)={}'.format(n))







plt.xlabel('Time Difference [ns]')

plt.ylabel('Coincidences')



plt.plot(x_data,y_data_bg_corrected,color='black',label='Simulation',alpha=0.6)

plt.plot(x_,gaussian(x_,popt[0],popt[1],popt[2]),label='Fitted Gaussian',color='red',linestyle='dotted')

plt.xlim(-7.5e-9,7.5e-9)

plt.legend()

plt.savefig('histogram.jpg')

plt.show()


