import pymol
import os
import pickle
import subprocess
from pymol import cgo, cmd

from .calcPorosity import *

def angleToColour(angle, minAngle, maxAngle):
    """returns normalised values for a spectrum of blue to red"""
    norm = (angle - minAngle) / (maxAngle - minAngle)
    norm = max(0, min(1, norm))
    r = norm
    g = 0
    b = 1 - norm
    return r, g, b

def coordToColour(coordNum, minCoord = 0, maxCoord = 2):
    """returns normalised values for a spectrum of"""
    norm = (coordNum - minCoord) / (maxCoord - minCoord)
    norm = max(0, min(1, norm))
    r = 0
    g = norm
    b = 1- norm
    return r, g, b

def loadForceVis(barrier, refoldBarrier, runNum ,vf, numMol, boxLength, particles, unfoldedMols, angles, bondsPerAtom, prob, suffix, simType, strainFraction):
    """visualises trajectories over time according to unfolding degree of molecule."""
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    pymol.finish_launching()
    # runNum = 10
    parRadius = 2 ** (1/6) * 0.4
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    #pores = calcPorosity(barrier, refoldBarrier, runNum, vf, numMol, boxLength)
    frames = len(timesteps)
    minAngle = min(min(angles.values()))
    maxAngle = max(max(angles.values()))
    print(minAngle, maxAngle)
    for frame in range(frames):
        frame_cgo = []
        
        for num, p in enumerate(particles):
            x = particles[num].properties[frame, 1]
            y = particles[num].properties[frame, 2]
            z = particles[num].properties[frame, 3]
            molID = int(particles[num].properties[frame, 5])
            
            if molID < len(angles[frame]):
                angle = angles[frame][molID]
        
            r, g, b = angleToColour(angle, minAngle, maxAngle)
            
            # if 10 < z < 25:
            frame_cgo.extend([cgo.COLOR, r, g, b, cgo.SPHERE, x, y, z, parRadius])
        # for pore in pores:
        #     x, y, z = pore["pos"]
        #     poreRadius = pore["size"]

        #     r, g, b = 0, 0.5, 0
        #     frame_cgo.extend([cgo.COLOR, r, g, b, cgo.SPHERE, x, y, z, poreRadius])
        cmd.load_cgo(frame_cgo, f"angle_unfold{barrier}_Run{runNum}_Vf{vf}_mol{numMol}", state = frame + 1)

    cmd.zoom("all", buffer = 5)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    # framesDir = os.path.abspath("./frames")
    # saveMovie(vf, numMol, runNum, timesteps)
    # framesToMp4(f"Run{runNum}_Vf{vf}_mol{numMol}.mp4", framesDir, fps = 10)
    return

def loadCoordVis(barrier, refoldBarrier, runNum, vf, numMol, particles, timesteps, parCoordination):
    """visualises trajectories in pymol and colours according to extent of molecule coordination."""
    pymol.finish_launching()
    
    parRadius = 2 ** (1/6) * 0.4
    
    for frame in range(len(timesteps)):
        frame_cgo = []
        
        for num, p in enumerate(particles):
            x = particles[num].properties[frame, 1]
            y = particles[num].properties[frame, 2]
            z = particles[num].properties[frame, 3]
            molID = int(particles[num].properties[frame, 5])
            
            if molID < len(parCoordination[frame]):
                coord = parCoordination[frame][molID]
        
            r, g, b = coordToColour(coord)
            
            frame_cgo.extend([cgo.COLOR, r, g, b, cgo.SPHERE, x, y, z, parRadius])

        cmd.load_cgo(frame_cgo, f"coordination_unfold{barrier}_Run{runNum}_Vf{vf}_mol{numMol}", state = frame + 1)
        
    cmd.zoom("all", buffer = 5)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    # framesDir = os.path.abspath("./frames")
    # saveMovie(vf, numMol, runNum, timesteps)
    # framesToMp4(f"Run{runNum}_Vf{vf}_mol{numMol}.mp4", framesDir, fps = 10)
    return
