import json,sys,os

filename=sys.argv[1]
with open(filename) as f: d=json.load(f)["active"]["MGxss"]

mats=[x[0] for x in d["MaterialBins"]]
ng=len(d["EnergyBounds"][0])
out=os.path.splitext(filename)[0]+"_xs"
os.makedirs(out,exist_ok=True)

def mean(x): return [v[0] for v in x]

for m in range(len(mats)):
    capture=mean(d["capture"][m])
    fission=mean(d["fission"][m])
    outscatter=mean(d["transportOutScatter"][m])
    total=[a+b+c for a,b,c in zip(capture,fission,outscatter)]
    nu=mean(d["nu"][m])
    chi=mean(d["chi"][m])
    p0=mean(d["P0"][m])
    p1=mean(d["P1"][m])

    s=f"// These XSs were generated from {os.path.basename(filename)}\n\n"
    s+=f"numberOfGroups {ng};\n"
    for name,x in [("total",total),("capture",capture),("fission",fission),("nu",nu),("chi",chi)]:
        s+=f"{name} ( "+" ".join(f"{v:.7E}" for v in x)+");\n"

    s+="\nscatteringMultiplicity (\n"
    for i in range(ng):
        s+=" ".join("1.0" for _ in range(ng))+"\n"
    s+=");\n\n"

    s+="P0 (\n"
    for i in range(ng):
        s+=" ".join(f"{v:.7E}" for v in p0[i*ng:(i+1)*ng])+"\n"
    s+=");\n\n"

    s+="P1 (\n"
    for i in range(ng):
        s+=" ".join(f"{v:.7E}" for v in p1[i*ng:(i+1)*ng])+"\n"
    s+=");\n"

    with open(os.path.join(out,mats[m]),"w") as f: f.write(s)
