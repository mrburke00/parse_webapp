# parse_webapp

# Setup (In Progress)
Clone Repo
```
cd parse_webapp/
```
Install vritualenv
```
pip install virtualenv
```
Find python path (i.e /usr/local/bin/python3.6) and start virtualenv 
```
virtualenv -p /usr/local/bin/python3.6 cmd2web_env
source cmd2web_env/bin/activate
```
Run configure script
```
bash parse_src/configure.sh 
```
Change `command` path in `cmd2web/ex_configs/tabix_config.json
```
        "name" : "parse",
        ...
        "command":
        [
            "gfortran -o /home/user/parse_webapp/cmd2web/src/web_src/static/js/Parse.exe /home/user/parse_webapp/cmd2web/src/web_src/static/js/Parse.f",
            "&& /home/user/parse_webapp/cmd2web/src/web_src/static/js/./Parse.exe ",
            "$sequence"
        ],
```
Start Server (need to change a few paths in tabix_config.json)
```
python cmd2web/src/server.py --config cmd2web/ex_configs/tabix_config.json 
```
Navigate to http://127.0.0.1:8080/parse
Enter query sequence into search bar (note: hit enter to submit, button is broken right now)
