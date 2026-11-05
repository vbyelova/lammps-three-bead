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


def loadBondForceVis(barrier, refoldBarrier, runNum ,vf, numMol, boxLength, particles, bondsPerAtom, prob, suffix, simType, strainFraction):
    """loads in bond information and visualises them as cylinders"""
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    pymol.finish_launching()
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
    frames = len(timesteps)
    with open(f"../runs/{conditions}/{filename}/{prefix}analysis/bondInfo.pkl", "rb") as f:
        bondInfo = pickle.load(f)
    allBonds = np.vstack(list(bondInfo.values()))
    mags = np.linalg.norm(allBonds[:, 3:6], axis=1)
    minForce, maxForce = mags.min(), mags.max()
    for frame in range(frames):
        frame_cgo = []
        bonds = bondInfo[frame]
        for b in range(len(bonds)):
            bond = bonds[b]
            atom1index = int(bond[1]) - 1
            atom2index = int(bond[2]) - 1
            largestBondForce = np.linalg.norm(bond[3:6])
            colour = angleToColour(largestBondForce, minForce, maxForce)

            x1 = particles[atom1index].properties[frame, 1]
            y1 = particles[atom1index].properties[frame, 2]
            z1 = particles[atom1index].properties[frame, 3]

            x2 = particles[atom2index].properties[frame, 1]
            y2 = particles[atom2index].properties[frame, 2]
            z2 = particles[atom2index].properties[frame, 3]            
            # if 10 < z < 25:
            frame_cgo.extend([cgo.CYLINDER, x1, y1, z1, x2, y2, z2, 0.5, *colour, *colour])
        # for pore in pores:
        #     x, y, z = pore["pos"]
        #     poreRadius = pore["size"]

        #     r, g, b = 0, 0.5, 0
        #     frame_cgo.extend([cgo.COLOR, r, g, b, cgo.SPHERE, x, y, z, poreRadius])
        cmd.load_cgo(frame_cgo, f"angle_unfold{barrier}_Run{runNum}_Vf{vf}_mol{numMol}", state = frame + 1)

    cmd.zoom("all", buffer = 5)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
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
