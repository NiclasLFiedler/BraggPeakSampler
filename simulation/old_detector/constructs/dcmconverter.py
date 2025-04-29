import pydicom
import numpy as np

def get_hu_values(dicom_file):
    dicom_data = pydicom.dcmread(dicom_file)
    image = dicom_data.pixel_array
    intercept = dicom_data.RescaleIntercept
    slope = dicom_data.RescaleSlope
    hu_values = image * slope + intercept
    return hu_values

filenames_source = [f"dicom_source/DICOM_Waterphantom{i}.dcm" for i in range(1, 782)]
filenames_target = [f"dicom_conv/DICOM_Waterphantom{i}.txt" for i in range(1, 782)]

dicom_data = pydicom.dcmread(filenames_source[0])
print(f"{dicom_data.PixelSpacing}")

for index,file in enumerate(filenames_source):
    hu_values = get_hu_values(file)
    with open(filenames_target[index], 'w') as file:
        # Step 3: Write each row of the matrix to the file
        for row in hu_values:
            # Convert each element to string and join them with spaces
            file.write(' '.join(map(str, row)) + '\n')