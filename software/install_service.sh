set -xe
cp dmbdi.service /etc/systemd/system/
systemctl daemon-reload
systemctl start dmbdi.service
systemctl enable dmbdi.service
