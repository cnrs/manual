github中文社区项目：  
https://www.githubs.cn/awesome  


samba：

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
        
[wangk]
        # wangk 用户的共享的目录
        path=/samba/wangk
        # 禁止匿名访问
        public=no
        # 是否可写
        writable=yes
        # 目录可写的用户组
        write list=@user
        # 访问目录的用户
        valid users=wangk

useradd wangk
smbpasswd -a wangk  #然后输入两次密码就可


mkdir -p /data/share
cd /data
chmod -R 775 share
chown -R wangk:user share


[root@cs home]# testparm
Load smb config files from /etc/samba/smb.conf
Loaded services file OK.
Server role: ROLE_STANDALONE

Press enter to see a dump of your service definitions

samba 服务器启动会后，默认会监听 139 和 445 端口，可以通过下面的命令查看 samba 服务器侦听的端口
[root@cs home]# netstat -an4p | grep smbd | grep LISTEN
tcp     0    0 0.0.0.0:139     0.0.0.0:*      LISTEN      23370/smbd          
tcp     0    0 0.0.0.0:445     0.0.0.0:*      LISTEN      23370/smbd

如果 samba 服务器所在的机器上开启了防火墙服务，则需要开放 139 和 445 端口，然后重启防火墙服务
[root@cs ~]# firewall-cmd --zone=public --add-port=139/tcp --permanent
success
[root@cs ~]# firewall-cmd --zone=public --add-port=445/tcp --permanent
success
[root@cs ~]# systemctl restart firewalld

service smb status
systemctl enable smb
systemctl start smb
```
