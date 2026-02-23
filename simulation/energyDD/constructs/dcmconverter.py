import pydicom
import numpy as np
import os

def get_hu_values(dicom_file):
    dicom_data = pydicom.dcmread(dicom_file)
    image = dicom_data.pixel_array
    intercept = dicom_data.RescaleIntercept
    slope = dicom_data.RescaleSlope
    hu_values = image * slope + intercept
    return hu_values

def convert(fpmod, fthickness):
    folder_path = f"dicom_source/{fpmod}mu_{fthickness}mm"
    num_files = len([f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))])
    filenames_source = [f"dicom_source/{fpmod}mu_{fthickness}mm/DICOM_Waterphantom{i}.dcm" for i in range(1, num_files+1)]
    filenames_target = [f"dicom_conv/{fpmod}mu_{fthickness}mm/DICOM_Waterphantom{i}.txt" for i in range(1, num_files+1)]

    dicom_data = pydicom.dcmread(filenames_source[0])
    print(f"{dicom_data.PixelSpacing}")

    for index,file in enumerate(filenames_source):
        hu_values = get_hu_values(file)
        with open(filenames_target[index], 'w') as file:
            # Step 3: Write each row of the matrix to the file
            for row in hu_values:
                # Convert each element to string and join them with spaces
                file.write(' '.join(map(str, row)) + '\n')

pmod = {100,200,300,400,500,600,700,800}
thickness = {50,100,150,200}

for p in pmod:
    for t in thickness:
        convert(p, t)
