import pymol
import os
import subprocess
from pymol import cgo, cmd

def saveMovie(vf, numMol, runNum, timesteps):
    cmd.set("ray_trace_frames", 0)
    cmd.set("movie_fps", 10)
    cmd.refresh()
    os.makedirs("./frames", exist_ok = True)
    cmd.mpng("./frames/frame_")
    
    import time
    time.sleep(2)
    return
    
def framesToMp4(output, framesDir, fps = 10):
    inputPattern = os.path.join(framesDir, "frame_%04d.png")
    result = subprocess.run(["ffmpeg", "-framerate", str(fps), "-i", "./frames/frame_%04d.png",
                            "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", output],
                            capture_output = True, text = True)
    print(result.stdout)
    print(result.stderr)
    return

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

def loadForceVis(barrier, refoldBarrier, runNum ,vf, numMol, particles, timesteps, unfoldedMols, angles):
    pymol.finish_launching()
    
    parRadius = 2 ** (1/6) * 0.4
    
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
            
            frame_cgo.extend([cgo.COLOR, r, g, b, cgo.SPHERE, x, y, z, parRadius])

        cmd.load_cgo(frame_cgo, f"angle_unfold{barrier}_Run{runNum}_Vf{vf}_mol{numMol}", state = frame + 1)
        
    cmd.zoom("all", buffer = 5)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    # framesDir = os.path.abspath("./frames")
    # saveMovie(vf, numMol, runNum, timesteps)
    # framesToMp4(f"Run{runNum}_Vf{vf}_mol{numMol}.mp4", framesDir, fps = 10)
    return

def loadCoordVis(barrier, refoldBarrier, runNum, vf, numMol, particles, timesteps, parCoordination):

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
