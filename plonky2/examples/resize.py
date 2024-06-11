import math
import numpy as np
import json
import sys

def bilinear_circom(image, hOrig, wOrig, hNew, wNew):

  resized = np.empty([hNew, wNew, 3])
  positiveRemainder = np.empty([hNew, wNew, 3])
  negativeRemainder = np.empty([hNew, wNew, 3])

  for i in range(hNew):
    for j in range(wNew):
      for k in range(3):
           x_l2 = (int)((wOrig - 1) * j / (wNew - 1))
           y_l2 = (int)((hOrig - 1) * i / (hNew - 1))
           x_h2 = x_l2 if x_l2 * (wNew - 1) == (wOrig - 1) * j else x_l2 + 1
           y_h2 = y_l2 if y_l2 * (hNew - 1) == (hOrig - 1) * i else y_l2 + 1

           print(x_l2, y_l2, x_h2, y_h2)

           xRatioWeighted = ((wOrig - 1) * j) - (wNew - 1)*x_l2
           yRatioWeighted = ((hOrig - 1) * i) - (hNew - 1)*y_l2

           a = image[y_l2, x_l2, k]
           b = image[y_l2, x_h2, k]
           c = image[y_h2, x_l2, k]
           d = image[y_h2, x_h2, k]

           s = a * (wNew - 1 - xRatioWeighted) * (hNew - 1 - yRatioWeighted) \
                  + b * xRatioWeighted * (hNew - 1 - yRatioWeighted) \
                  + c * yRatioWeighted * (wNew - 1 - xRatioWeighted) \
                  + d * xRatioWeighted * yRatioWeighted

           new = (round)( s / ((wNew - 1) * (hNew - 1)))

           r = (round(new * (wNew - 1) * (hNew - 1) - s))
           
           resized[i][j][k] = new
           if r > 0:
                positiveRemainder[i][j][k] = r
                negativeRemainder[i][j][k] = 0
           else:
                negativeRemainder[i][j][k] = abs(r)
                positiveRemainder[i][j][k] = 0
 
hOrig = int(sys.argv[1])
wOrig = int(sys.argv[2])
hNew = int(sys.argv[3])
wNew = int(sys.argv[4])
x = np.random.randint(256, size=(hOrig, wOrig, 3))
#print(x)
bilinear_circom(x, hOrig, wOrig, hNew, wNew)