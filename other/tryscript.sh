vboxmanage snapshot omegues restore ooo
vboxmanage startvm omegues --type headless
sleep 5
ssh omegues uname -sr
vboxmanage controlvm omegues poweroff
