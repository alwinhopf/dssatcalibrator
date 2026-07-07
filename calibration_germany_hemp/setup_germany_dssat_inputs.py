"""Install the derived Germany hemp inputs into the local DSSAT48Hemp tree.

The PDF extraction script writes derived weather, soil, observations, and
management tables under ``derived/``. This setup step creates the runnable DSSAT
side: two FileX files, provisional cultivar/ecotype rows, the central weather
file, and the reported-soil profile.
"""

from __future__ import annotations

from datetime import date, timedelta
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DERIVED = HERE / "derived"
DSSAT = ROOT / "DSSAT48Hemp"
HEMP_DIR = DSSAT / "Hemp"
GENO_DIR = DSSAT / "Genotype"
WEATHER_DIR = DSSAT / "Weather"
SOIL_DIR = DSSAT / "Soil"

SITE_ID = "GEMQ2018"
HARVEST = date(2018, 9, 19)
IRRIGATION = date(2018, 5, 30)

EXPERIMENTS = [
    {
        "exp_id": "GEMQ18S1",
        "title": "Santhica 27",
        "cultivar": "DE0001",
        "cname": "Santhica27",
        "pdate": date(2018, 5, 4),
        "ppop": 200.0,
        "ppoe": 149.0,
        "area_ha": 0.12,
        "anchor_cul": "IB0007",
        "anchor_eco": "HM0005",
        "new_eco": "HMDE01",
        "cul_name": "Santhica27Scaffold",
        "eco_name": "Santhica27 Fiber",
    },
    {
        "exp_id": "GEMQ18I1",
        "title": "Ivory",
        "cultivar": "DE0002",
        "cname": "Ivory",
        "pdate": date(2018, 5, 22),
        "ppop": 200.0,
        "ppoe": 60.0,
        "area_ha": 0.08,
        "anchor_cul": "IB0005",
        "anchor_eco": "HM0001",
        "new_eco": "HMDE02",
        "cul_name": "IvoryScaffold",
        "eco_name": "Ivory Fiber",
    },
]


def yyddd(d: date) -> int:
    return int(f"{d.year % 100:02d}{d.timetuple().tm_yday:03d}")


def _line_for_code(lines: list[str], code: str) -> str:
    try:
        return next(ln for ln in lines if ln.startswith(code) and not ln.lstrip().startswith("!"))
    except StopIteration as exc:
        raise RuntimeError(f"Could not find required DSSAT row {code!r}") from exc


def _replace_slice(line: str, start: int, end: int, value: str) -> str:
    width = end - start
    return line[:start] + value[:width].ljust(width) + line[end:]


def install_genotypes() -> None:
    cul_path = GENO_DIR / "HMGRO048.CUL"
    eco_path = GENO_DIR / "HMGRO048.ECO"
    cul_lines = cul_path.read_text(errors="replace").splitlines()
    eco_lines = eco_path.read_text(errors="replace").splitlines()

    added_cul: list[str] = []
    added_eco: list[str] = []
    for exp in EXPERIMENTS:
        if not any(ln.startswith(exp["cultivar"]) for ln in cul_lines):
            row = _line_for_code(cul_lines, exp["anchor_cul"])
            row = _replace_slice(row, 0, 6, exp["cultivar"])
            row = _replace_slice(row, 7, 25, exp["cul_name"])
            row = _replace_slice(row, 30, 36, exp["new_eco"])
            added_cul.append(row)
        if not any(ln.startswith(exp["new_eco"]) for ln in eco_lines):
            row = _line_for_code(eco_lines, exp["anchor_eco"])
            row = _replace_slice(row, 0, 6, exp["new_eco"])
            row = _replace_slice(row, 7, 24, exp["eco_name"])
            added_eco.append(row)

    if added_cul:
        cul_lines.extend(["", "!Germany Marquardt/Potsdam scaffolds; provisional before calibration."])
        cul_lines.extend(added_cul)
        cul_path.write_text("\n".join(cul_lines) + "\n", encoding="utf-8")
    if added_eco:
        eco_lines.extend(["", "!Germany Marquardt/Potsdam scaffolds; provisional before calibration."])
        eco_lines.extend(added_eco)
        eco_path.write_text("\n".join(eco_lines) + "\n", encoding="utf-8")


def install_weather() -> None:
    src = DERIVED / "weather_dwd" / f"{SITE_ID}.WTH"
    dst = WEATHER_DIR / f"{SITE_ID}.WTH"
    lines = src.read_text(encoding="utf-8").splitlines()
    for i, line in enumerate(lines):
        if line.startswith("$WEATHER DATA"):
            lines[i] = "$WEATHER DATA: DWD Potsdam station 03987 for Germany hemp site"
        if i > 0 and lines[i - 1].lstrip().startswith("@ INSI"):
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                lines[i] = f"  GEMQ  {parts[1]}"
            break
    dst.write_text("\n".join(lines) + "\n", encoding="utf-8")


def install_soil() -> None:
    src = DERIVED / "soil" / f"{SITE_ID}_reported_approx.SOL"
    soil_sol = SOIL_DIR / "SOIL.SOL"
    profile = src.read_text(encoding="utf-8").strip()
    current = soil_sol.read_text(errors="replace") if soil_sol.exists() else ""
    if f"*{SITE_ID}" not in current:
        with soil_sol.open("a", encoding="utf-8") as fh:
            if current and not current.endswith("\n"):
                fh.write("\n")
            fh.write("\n")
            fh.write(profile)
            fh.write("\n")


def hmx_text(exp: dict) -> str:
    pdate = yyddd(exp["pdate"])
    fdate = yyddd(exp["pdate"] - timedelta(days=1))
    idate = yyddd(IRRIGATION)
    hdate = yyddd(HARVEST)
    pfrst = yyddd(exp["pdate"] - timedelta(days=5))
    plast = yyddd(exp["pdate"] + timedelta(days=10))
    name = exp["title"]
    icname = "75% relative soil water, 100 kg mineral N/ha"

    return f"""*EXP.DETAILS: {exp['exp_id']}HM Germany hemp {name}, Marquardt/Potsdam 2018

*GENERAL
@PEOPLE
Gabriele Gusovius et al.; digitized by dssatcalibrator scaffold
@ADDRESS
Leibniz Institute for Agricultural Engineering and Bioeconomy, Potsdam, Germany
@SITE
Marquardt/Potsdam field site, Germany
@ PAREA  PRNO  PLEN  PLDR  PLSP  PLAY HAREA  HRNO  HLEN  HARM.........
    -99   -99   -99   -99   -99   -99   -99   -99   -99   -99

*TREATMENTS                        -------------FACTOR LEVELS------------
@N R O C TNAME.................... CU FL SA IC MP MI MF MR MC MT ME MH SM
 1 1 1 0 {name[:25]:<25}  1  1  0  1  1  1  1  0  0  0  0  1  1

*CULTIVARS
@C CR INGENO CNAME
 1 HM {exp['cultivar']} {exp['cname']}

*FIELDS
@L ID_FIELD WSTA....  FLSA  FLOB  FLDT  FLDD  FLDS  FLST SLTX  SLDP  ID_SOIL    FLNAME
 1 {SITE_ID} {SITE_ID}   -99     0 IB000     0     0 00000 -99    200  {SITE_ID}   -99
@L ...........XCRD ...........YCRD .....ELEV .............AREA .SLEN .FLWR .SLAS FLHST FHDUR
 1          12.9608        52.4672        81              {exp['area_ha']:.2f}   -99   -99   -99   -99   -99

*SOIL ANALYSIS
@A SADAT  SMHB  SMPX  SMKE  SANAME
 1 {fdate}   -99   -99   -99  Reported low loamy sand SL2
@A  SABL  SADM  SAOC  SANI SAPHW SAPHB  SAPX  SAKE  SASC
 1   200   -99   .70   -99   -99   -99   -99   -99   -99

*INITIAL CONDITIONS
@C   PCR ICDAT  ICRT  ICND  ICRN  ICRE  ICWD ICRES ICREN ICREP ICRIP ICRID ICNAME
 1    HM {fdate}     0   -99     1     1   -99     0     0   -99     0    10 {icname}
@C  ICBL  SH2O  SNH4  SNO3
 1     5  .219   0.5   2.9
 1    15  .219   0.5   2.9
 1    30  .219   0.5   2.9
 1    60  .219   0.5   2.9
 1   100  .219   0.5   2.9
 1   200  .219   0.5   2.9

*PLANTING DETAILS
@P PDATE EDATE  PPOP  PPOE  PLME  PLDS  PLRS  PLRD  PLDP  PLWT  PAGE  PENV  PLPH  SPRL                        PLNAME
 1 {pdate}   -99 {exp['ppop']:5.1f} {exp['ppoe']:5.1f}     S     R  12.5     0   2.0   -99   -99   -99   -99   -99                        200 seeds/m2, observed emergence stand

*IRRIGATION AND WATER MANAGEMENT
@I  EFIR  IDEP  ITHR  IEPT  IOFF  IAME  IAMT IRNAME
 1  1.00   -99   -99   -99   -99   -99   -99 Sprinkler
@I IDATE  IROP IRVAL
 1 {idate} IR004    10

*FERTILIZERS (INORGANIC)
@F FDATE  FMCD  FACD  FDEP  FAMN  FAMP  FAMK  FAMC  FAMO  FOCD FERNAME
 1 {fdate} FE001 AP002     1    70     0     0   -99   -99   -99 Calcium ammonium nitrate, assumed day before planting

*ENVIRONMENT MODIFICATIONS
@E ODATE EDAY  ERAD  EMAX  EMIN  ERAIN ECO2  EDEW  EWIND ENVNAME
 1 {pdate} A   0 A   0 A   0 A   0 A 0.0 A   0 A   0 A   0

*HARVEST DETAILS
@H HDATE  HSTG  HCOM HSIZE   HPC  HBPC HNAME
 1 {hdate} GS000   -99     A   -99   -99 Water paper oven-dry harvest date

*SIMULATION CONTROLS
@N GENERAL     NYERS NREPS START SDATE RSEED SNAME.................... SMODEL
 1 GE              1     1     S {fdate}  2150 Germany 2018 hemp        CRGRO
@N OPTIONS     WATER NITRO SYMBI PHOSP POTAS DISES  CHEM  TILL   CO2
 1 OP              Y     Y     N     N     N     N     N     N     M
@N METHODS     WTHER INCON LIGHT EVAPO INFIL PHOTO HYDRO NSWIT MESOM MESEV MESOL
 1 ME              M     M     E     R     S     L     R     1     P     R     2
@N MANAGEMENT  PLANT IRRIG FERTI RESID HARVS
 1 MA              R     R     R     N     R
@N OUTPUTS     FNAME OVVEW SUMRY FROPT GROUT CAOUT WAOUT NIOUT MIOUT DIOUT VBOSE CHOUT OPOUT FMOPT
 1 OU              N     Y     Y     1     Y     Y     Y     Y     N     N     Y     N     N     A

@  AUTOMATIC MANAGEMENT
@N PLANTING    PFRST PLAST PH2OL PH2OU PH2OD PSTMX PSTMN
 1 PL          {pfrst} {plast}    40   100    30    40    10
@N IRRIGATION  IMDEP ITHRL ITHRU IROFF IMETH IRAMT IREFF
 1 IR             30    70   100 IB001 IB001    10   .75
@N NITROGEN    NMDEP NMTHR NAMNT NCODE NAOFF
 1 NI             30    50    25 IB001 IB001
@N RESIDUES    RIPCN RTIME RIDEP
 1 RE            100     1    20
@N HARVEST     HFRST HLAST HPCNP HPCNR
 1 HA              0 18365   100     0

"""


def install_filex() -> None:
    HEMP_DIR.mkdir(parents=True, exist_ok=True)
    for exp in EXPERIMENTS:
        path = HEMP_DIR / f"{exp['exp_id']}.HMX"
        path.write_text(hmx_text(exp), encoding="utf-8")


def copy_observation_csv_to_hemp_dir() -> None:
    """Keep a convenience copy next to FileX without making DSSAT depend on it."""
    src = DERIVED / "dssatcalibrator_observations_long.csv"
    if src.exists():
        shutil.copy(src, HEMP_DIR / "GEMQ18_observations_long.csv")


def main() -> None:
    install_genotypes()
    install_weather()
    install_soil()
    install_filex()
    copy_observation_csv_to_hemp_dir()
    print("Installed Germany hemp DSSAT inputs:")
    print(f"  FileX: {HEMP_DIR / 'GEMQ18S1.HMX'}")
    print(f"         {HEMP_DIR / 'GEMQ18I1.HMX'}")
    print(f"  Weather: {WEATHER_DIR / (SITE_ID + '.WTH')}")
    print(f"  Soil profile: {SOIL_DIR / 'SOIL.SOL'} (*{SITE_ID})")
    print(f"  Genotypes: DE0001/HMDE01 and DE0002/HMDE02 if they were missing")


if __name__ == "__main__":
    main()
