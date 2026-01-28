import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import MDAnalysis as mda
import numpy as np
from pathlib import Path
import ctleelab_plothelper.plothelpers as ph
from PIL import Image

configuration = "strain2"
movietype = "lipid-type"
views = ["side","top"]
lipid_compositions = {
    "DOPC": [2, 0],
    "DPPC": [2, 1],
    "DOPC-DOPA": [1, 0],
    "DPPC-DPPA": [1, 1],
    "DOPC-DOPS": [0, 0],
    "DPPC-DPPS": [0, 1],
}


# setting paths
script_dir = Path(__file__).resolve().parent
systems_path = Path(script_dir.parent / "Figures"/"MD-images"/movietype/configuration)
compiled_path = systems_path/"compiled"/"all-systems"
compiled_path.mkdir(exist_ok=True)


for view in views: 
    for number in range(30): #loops through frames
        image_number = f"{number:04d}"
        with plt.style.context(["ctleelab_plothelper.base", "ctleelab_plothelper.light"]):
            fig, ax = ph.fixed_size_subplots(3,2,subwidth=3,subheight=2.5
            ,wmargin=0.05, hmargin=0.002, rmargin_scale = 0.1,tmargin_scale = 0.2)
            for composition, position in lipid_compositions.items():
                folder_system = systems_path/composition/view
                file = folder_system/f"frame{image_number}.png"

                # Read PNG as RGBA
                img = Image.open(file).convert("RGBA")
                data = np.array(img)

                white = (data[:, :, 0:3] == 255).all(axis=2) #determines which pixels are white
                data[white, 3] = 0  # alpha = 0 for white
                alpha = data[:, :, 3]

                #remove the transparent background 
                rows = np.any(alpha > 0, axis=1)
                cols = np.any(alpha > 0, axis=0)
                top, bottom = np.where(rows)[0][[0, -1]]
                left, right = np.where(cols)[0][[0, -1]]

                cropped_img = img.crop((left, top, right+1, bottom+1))


                i, j = position
                ax[i,j].imshow(cropped_img)
                ax[i, j].axis("off")   
                ax[i,j].set_title(composition, fontsize = 16)
            plt.savefig(compiled_path/f"{view}-frame{number}.png", transparent=True)
        plt.close(fig)
