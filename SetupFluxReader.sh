
source /cvmfs/nova.opensciencegrid.org/externals/setup
setup cmake v3_18_2


if(env | grep -q ^DK2NU)
then
  echo "Dk2Nu has already been setup."
  echo "If the version is not correct, this may cause errors."
  echo "If this occurs, try starting a fresh terminal session."
else
  setup dk2nudata v01_08_00g -q debug:e19
fi

if(env | grep -q ^GENIEXSECPATH)
then
  echo "genie_xsec has already been setup."
  echo "If the version is not correct, this may cause errors."
  echo "If this occurs, try starting a fresh terminal session."
else
  setup genie_xsec v2_10_6 -q defaultplusccmec
fi

export FLUXREADER_PRIV="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "FLUXREADER_PRIV set to ${FLUXREADER_PRIV}"
