import cv2
import numpy as np


import io, struct, time, sys

fname = 'C:/Users/Zahra/Documents/GitHub/Data/gabardine_od_supertile.vol'
flipped = True
prctle = 0.95
processed = False

with io.open(fname, mode='rb') as fin:
    assert fin.read(3) == 'VOL'
    assert fin.read(1) == '\x03'
    assert struct.unpack('I', fin.read(4)) == (1,)

    time.clock()
    sz = struct.unpack('3I', fin.read(12))
    ch = struct.unpack('I', fin.read(4))
    assert(len(ch) > 0 and ch[0] == 1)

    pMin = np.array(struct.unpack('3f', fin.read(12)))
    pMax = np.array(struct.unpack('3f', fin.read(12)))
    extent = pMax - pMin

    sys.stdout.write('Loading: ')
    sys.stdout.flush()
    A = np.empty((sz[2], sz[1], sz[0]))
    n = sz[1]*sz[0]
    for i in range(0, sz[2]):
        A[i, :, :] = np.reshape(np.asarray(struct.unpack('%df' % n, fin.read(4*n))), (sz[1], sz[0]), 'C')
        if i % 50 == 0:
            sys.stdout.write('\rLoading: %d/%d' %(i, sz[2]))
            sys.stdout.flush()
    sys.stdout.write('\rLoaded in %.2f secs.\n' % time.clock())

    A1 = np.zeros((sz[2], sz[1], sz[0], 3))
    A1[:, :, :, 0] = np.abs(A)
    A1[:, :, :, 1] = np.abs(A)
    A1[:, :, :, 2] = np.abs(A)

    if flipped:
        threshold = sz[2] - 1
    else:
        threshold = 0
    height = -1.0*np.ones((sz[1], sz[0]))

    mode = 1
    x, y, z = 0, 0, 0
    while True:
        if mode == 1:
            B = np.copy(A1[:, :, x, :])
            cv2.line(B, (0, threshold), (sz[1], threshold), (0, 255, 0))
            pts = []
            for i in range(0, sz[1]):
                if height[i, x] > 0.0:
                    pts.append([i, int(height[i, x])])
            pts = np.array(pts, np.int32).reshape((-1,1,2))
            if len(pts) > 1:
                cv2.polylines(B, [pts], False, (0,255,255))
        elif mode == 2:
            B = np.copy(A1[:, y, :, :])
            cv2.line(B, (0, threshold), (sz[0], threshold), (0, 255, 0))
            pts = []
            for i in range(0, sz[0]):
                if height[y, i] > 0.0:
                    pts.append([i, int(height[y, i])])
            pts = np.array(pts, np.int32).reshape((-1,1,2))
            if len(pts) > 1:
                cv2.polylines(B, [pts], False, (0,255,255))
        else:
            B = A1[z, :, :, :]
        cv2.imshow('Image', B)

        key = cv2.waitKey(0)
        if key == 27 or key < 0:
            break
        elif key == 2424832:
            if mode == 1:
                x = np.max([0, x - 1])
            elif mode == 2:
                y = np.max([0, y - 1])
            else:
                z = np.max([0, z - 1])
        elif key == 2555904:
            if mode == 1:
                x = np.min([sz[0] - 1, x + 1])
            elif mode == 2:
                y = np.min([sz[1] - 1, y + 1])
            else:
                z = np.min([sz[2] - 1, z + 1])
        elif key == 2490368:
            threshold = max(threshold - 1, 0)
        elif key == 2621440:
            threshold = min(threshold + 1, sz[2] - 1)
        elif key == ord('1'):
            mode = 1
        elif key == ord('2'):
            mode = 2
        elif key == ord('3'):
            mode = 3
        elif key == 13:
            # Process Volume
            sys.stdout.write('Processing ... ')
            sys.stdout.flush()

            height = -1.0*np.ones((sz[1], sz[0]))
            if flipped:
                for i in range(0, sz[1]):
                    for j in range(0, sz[0]):
                        v = np.cumsum(A[0 : threshold + 1, i, j])
                        if v[-1] > 0:
                            v /= v[-1]
                            k = threshold
                            while v[k] > prctle:
                                k -= 1
                            height[i, j] = k
            else:
                for i in range(0, sz[1]):
                    for j in range(0, sz[0]):
                        v = np.cumsum(A[threshold :, i, j])
                        if v[-1] > 0:
                            v /= v[-1]
                            k = threshold
                            while v[k] < 1.0 - prctle:
                                k += 1
                            height[i, j] = k

            sys.stdout.write('done\n')
            processed = True
        elif key == 32:
            # Save
            if processed:
                height_out = np.copy(height).astype(np.float32)
                height_out[height_out < 0] = np.min(height_out[height_out > 0])
                height_out = extent[2]*((height_out + 1.0)/sz[2])
                cv2.imwrite('heightmap.exr', height_out.astype(np.float32))

                with open('output.txt', 'w') as fout:
                    fout.write('<!-- Transform for the volume -->\n')
                    fout.write('<transform name="toWorld">\n')
                    fout.write('    <translate x="%.4f" y="%.4f" z="%.4f"/>\n' % (-0.5*(pMin[0] + pMax[0]), -0.5*(pMin[1] + pMax[1]), -pMin[2]))
                    fout.write('</transform>\n\n')

                    hybridZ = extent[2]*float(sz[2] - 1 - threshold)/sz[2]

                    fout.write('<!-- Transform for the bounding cube -->\n')
                    fout.write('<transform name="toWorld">\n')
                    fout.write('    <translate z="1"/>\n')
                    fout.write('    <scale x="%.4f" y="%.4f" z="%.4f"/>\n' % (0.5*extent[0], 0.5*extent[1], 0.5*hybridZ))
                    fout.write('    <translate z="%.4f"/>\n' % (extent[2]*(threshold + 1)/sz[2]))
                    fout.write('</transform>\n\n')

                    fout.write('<!-- Transform for the height map -->\n')
                    fout.write('<transform name="toWorld">\n')
                    fout.write('    <scale x="%.4f" y="%.4f"/>\n' % (extent[0]/2.0, extent[1]/2.0))
                    fout.write('</transform>\n\n')
            else:
                sys.stderr.write('Cannot save: volume has not been processed yet!\n')
        else:
            print key

sys.stdout.write('\rThe threshold is: %f .\n' % threshold)