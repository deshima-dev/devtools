SCP="/usr/bin/scp -Cp"
USER="nro"
HOST="nro.sakura.ne.jp"
LOG_DIR="/home/nro/tescam/cosmos3/aste/log/"
LOG_ID="20160728113318"

${SCP} ${USER}@${HOST}:${LOG_DIR}/obs/${LOG_ID}.obs .
${SCP} ${USER}@${HOST}:${LOG_DIR}/ant/${LOG_ID}.ant .
${SCP} ${USER}@${HOST}:${LOG_DIR}/weather/${LOG_ID}.wea .
