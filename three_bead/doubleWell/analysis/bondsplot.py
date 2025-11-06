import matplotlib.pyplot as plt
import pandas as pd

columns = ["time", "ssintra", "ssinter", "shinter", "hhinter"]
#time ssintra ssinter shinter hhinter
df = pd.read_csv("../output/nbonds.dat", skipinitialspace = True, delim_whitespace = True,  usecols = columns)
print(df.keys())

fig, axs = plt.subplots(2,2)

axs[0,0].plot(df.time, df.ssintra)
axs[0,0].set_title("intramol. sticker-sticker")
axs[0,1].plot(df.time, df.ssinter)
axs[0,1].set_title("intermol. sticker-sticker")
axs[1,0].plot(df.time, df.shinter)
axs[1,0].set_title("intermol. sticker-hinge")
axs[1,1].plot(df.time, df.hhinter)
axs[1,1].set_title("intermol. hinge-hinge")

for ax in axs.flat:
    ax.set(xlabel="time", ylabel="number of bonds")

#for ax in axs.flat:
#    ax.label_outer()

plt.show()
