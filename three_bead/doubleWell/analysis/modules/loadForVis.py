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
    
def framesToMp4(output, framesDir, fps = 10):
    inputPattern = os.path.join(framesDir, "frame_%04d.png")
    result = subprocess.run(["ffmpeg", "-framerate", str(fps), "-i", "./frames/frame_%04d.png",
                            "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", output],
                            capture_output = True, text = True)
    print(result.stdout)
    print(result.stderr)

def angleToColour(angle):
    norm = angle / 180
    r = norm
    g = 0
    b = 1 - norm
    return r, g, b

def loadForVis(barrier, refoldBarrier, runNum ,vf, numMol, particles, timesteps, unfoldedMols, angles):
    parRadius = 2 ** (1/6) * 0.3
    
    for frame in range(len(timesteps)):
        frame_cgo = []
        
        for num, p in enumerate(particles):
            x = particles[num].properties[frame, 1]
            y = particles[num].properties[frame, 2]
            z = particles[num].properties[frame, 3]
            molID = int(particles[num].properties[frame, 5]) - 1 
            
            if molID < len(angles[frame]):
                angle = angles[frame][molID]
        
            r, g, b = angleToColour(angle)
            
            frame_cgo.extend([cgo.COLOR, r, g, b, cgo.SPHERE, x, y, z, parRadius])
            
        cmd.load_cgo(frame_cgo, f"Run{runNum}_Vf{vf}_mol{numMol}", state = frame + 1)
        
    cmd.zoom("all", buffer = 5)
    framesDir = os.path.abspath("./frames")
    saveMovie(vf, numMol, runNum, timesteps)
    framesToMp4(f"Run{runNum}_Vf{vf}_mol{numMol}.mp4", framesDir, fps = 10)
    return