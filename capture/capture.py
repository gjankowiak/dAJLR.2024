import cv2
import numpy as np 
import matplotlib.pyplot as plt

import sys

###############################
#            CAPTURE          #
###############################

drawing=False # true if mouse is pressed
mode=True # if True, draw rectangle. Press 'm' to toggle to curve

dump_file = open("dump.txt", "w")

# mouse callback function
def interactive_drawing(event,x,y,flags,param):
    global ix,iy,drawing, mode

    if event==cv2.EVENT_LBUTTONDOWN:
        drawing=True
        ix,iy=x,y

    elif event==cv2.EVENT_MOUSEMOVE:
        if drawing==True:
            if mode==True:
                dump_file.write("%d %d\n" % (x, y))
                #print(x, y)
                cv2.circle(img,(x,y),1,(0,0,255),-1)
    elif event==cv2.EVENT_LBUTTONUP:
        drawing=False
        if mode==True:
            cv2.circle(img,(x,y),1,(0,0,255),-1)        


img = np.zeros((1024,1024,3), np.uint8)
cv2.namedWindow('Window')
cv2.setMouseCallback('Window', interactive_drawing)

print("Capturing to 'dump.txt', press ESC to finish", file=sys.stderr)

while(1):
    cv2.imshow('Window',img)
    k=cv2.waitKey(1)&0xFF
    if k==27:
        break

cv2.destroyAllWindows()
dump_file.close()

print("Done, you can now run postproc.py")
