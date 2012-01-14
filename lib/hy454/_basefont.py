

from os.path import exists

import numpy as np

from matplotlib.path import Path
import matplotlib.patches as patches

from freetype import *


class Basefont(object):

    def __init__(self, fontpath, charsize=48):
        if not exists(fontpath):
            raise ValueError('No valid font file specified')

        face = Face(fontpath)
        face.set_char_size(charsize * 64)
        self.face = face
        self.cache = {}

    def char_patch(self, char):
        # if we have the path in our cache, re-use it
        if char in self.cache:
            return patches.PathPatch(self.cache[char])

        self.face.load_char(char)
        slot = self.face.glyph
        outline = slot.outline

        start, end = 0, 0
        verts, codes = [], []

        for i in range(len(outline.contours)):
            end = outline.contours[i]
            points = outline.points[start:end+1]
            points.append(points[0])
            tags = outline.tags[start:end+1]
            tags.append(tags[0])

            verts.append(points[0])
            codes.append(Path.MOVETO)

            segments = [[points[0]]]
            for j in range(1, len(points)):
                segments[-1].append(points[j])
                if tags[j] & (1 << 0) and j < (len(points) - 1):
                    segments.append([points[j]])

            for segment in segments:
                if len(segment) == 2:
                    verts.extend(segment[1:])
                    codes.append(Path.LINETO)
                elif len(segment) == 3:
                    verts.extend(segments[1:])
                    codes.extend([Path.CURVE3, Path.CURVE3])
                else:
                    verts.append(segment[1])
                    codes.append(Path.CURVE3)
                    for i in range(1, len(segment)-2):
                        A, B = segment[i], segment[i + 1]
                        C = ((A[0] + B[0]) / 2., (A[1] + B[1]) / 2.)
                        verts.extend([C, B])
                        codes.extend([Path.CURVE3, Path.CURVE3])
                    verts.append(segment[-1])
                    codes.append(Path.CURVE3)

            start = end+1

        path = Path(verts, codes)
        self.cache[char] = path
        glyph = patches.PathPatch(path)

        return glyph
