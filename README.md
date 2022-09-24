samba：

https://www.cnblogs.com/12345huangchun/p/12268343.html

```
yum install samba


#/etc/samba/smb.conf

[global]
        workgroup = SAMBA           #设定 Samba Server 所要加入的工作组或者域。
        security = user             #设置用户访问Samba Server的验证方式，一共有四种验证方式
        passdb backend = tdbsam
        printing = cups
        printcap name = cups
        load printers = yes
        cups options = raw
[database]
        comment = share database     #这个是共享文件的描述
        path = /data/share           #设置共享文件夹的路径
        public = no                  #设置是否允许匿名访问
        writable = yes
        

useradd smb_root
smbpasswd -a smb_root  #然后输入两次密码就可

mkdir -p /data/share
cd /data
chmod -R 775 share
chown -R smb_root:smb_root share
```
