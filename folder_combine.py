import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def select_folder():
    root = tk.Tk()
    root.withdraw()
    folder_path = filedialog.askdirectory(
        title="Select Folder",
    )
    return folder_path

def load_csv_data(filepath, skip_header=3, skip_footer=4, usecols=None):
    return np.genfromtxt(filepath, delimiter=',', skip_header=skip_header, skip_footer=skip_footer, encoding='utf-8', usecols=usecols)

def combine_and_sum_frames(file_paths):
    combined_image = None
    filenames = []

    for path in file_paths:
        fname = os.path.splitext(os.path.basename(path))[0]
        filenames.append(fname)

        data = load_csv_data(path, usecols=range(1, 641))  # assumes 640 columns of data
        if data.ndim == 1:
            data = data[np.newaxis, :]  # ensure it's 2D if only one row

        if combined_image is None:
            combined_image = np.zeros_like(data, dtype=np.float64)

        combined_image += data.astype(np.float64)

    combined_name = "_".join(filenames)
    return combined_image, combined_name

def save_combined_image(image, combined_name, folder="combined_images"):
    os.makedirs(folder, exist_ok=True)
    save_path = os.path.join(folder, combined_name + ".png")
    plt.imsave(save_path, image, cmap='gray')
    print(f"Combined image saved to: {os.path.abspath(save_path)}")

# Main workflow
folder_path = select_folder()
if not folder_path:
    print("No folder selected, exiting.")
    exit()

file_paths = (glob.glob(os.path.join(folder_path, "*.csv")))
if not file_paths:
    print("No CSV files found in folder.")
    exit()

combined_image, combined_name = combine_and_sum_frames(file_paths)

plt.imshow(combined_image, cmap='gray')
plt.show()

save_combined_image(combined_image, combined_name)
