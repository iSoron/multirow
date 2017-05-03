import sys
import math

M=3
ZERO=1e-10

def angle_cmp(v):
    return math.atan2(v[0], v[1])

for fname in sys.argv[1:]:
    load(fname)
    nrays = len(pi)

    frat = f.apply_map(lambda x : x.nearby_rational(max_denominator=1000))
    vertices = [(f+R[i]/pi[i]).apply_map(lambda x: x.nearby_rational(max_denominator=1000)) for i in range(nrays)]
    poly = Polyhedron(vertices=vertices)

    print "Writing %s.tex..." % fname
    texfile = open("%s.tex" % fname, "w")
    texfile.write("\\draw[fill,fill opacity=0.1] ")
    for v in sorted(poly.vertices(), key=angle_cmp):
        texfile.write("(%s, %s) -- " % (v[0], v[1]))
    texfile.write("cycle;\n")

    #texfile.write("\\draw (%s, %s) circle [radius=1.5pt];\n" % (frat[0], frat[1]))

    for v in vertices:
        texfile.write("\\draw[dashed] (%s,%s) -- (%s,%s);\n" % (frat[0], frat[1], v[0], v[1]))

    texfile.close()

    print "Writing %s.png..." % fname
    p0 = sum([arrow(f, f + R[i] / pi[i], color='blue') for i in range(nrays)])
    pf = list_plot([f], color='red', pointsize=30, zorder=-1)
    p1 = poly.plot(color='blue', alpha=0.2)
    p3 = list_plot(
        [(x,y) for x in range(-M,M+1) for y in range(-M,M+1)],
        color='black', figsize=20, pointsize=30, zorder=1)

    output = "%s.png" % fname
    save(p0+pf+p3+p1, output, xmin=-M, ymin=-M, xmax=M, ymax=M)
