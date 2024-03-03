import json
from flask import Flask
from logging import debug
import logging
from flask import g, request, jsonify, session, send_file
from geqdsk import read
from mapping import mapping_core
from gensp import python_write
from utils import np_to_list
from flask_cors import CORS
from flask_socketio import SocketIO, emit
from utils import FigType
import time

import tempfile

app = Flask(__name__)
app.config['SECRET_KEY'] = "xsd9jgdyj<>"
CORS(app,supports_credentials=True,cors_allowed_origins="*",resources={r'/*': {'origins': '*'}})
socketio = SocketIO(app)
socketio.init_app(app, cors_allowed_origins="*")
xlogger = logging.getLogger('xmap')
logging.basicConfig(level=logging.WARNING)
flasklog = logging.getLogger('werkzeug')
flasklog.setLevel(logging.ERROR)


@socketio.on('connect')
def handle_test():
    crt = time.localtime()
    strtime = time.strftime("%Y-%m-%d %H:%M:%S",crt)
    print(f"[{strtime}] Connected!!!")
    emit('echo')

@socketio.on_error()
def error_handler(e):
    emit('flask_error', e)

@socketio.on('gfile')
def handle_gfile(params):
    file_str = params["data"]
    lsp = params["lsp"]
    lst = params["lst"]
    psimax_ratio = 0.995
    file_obj = tempfile.TemporaryFile(mode="w+")
    file_obj.write(file_str)
    file_obj.flush()
    file_obj.seek(0)
    try:
        gfile_data = read(file_obj)
        map_data, check_data = mapping_core(gfile_data, lsp, lst, psimax_ratio, figs=FigType.web, nR=200)
        new_check_data = np_to_list(check_data)
    except Exception as e:
        crt = time.localtime()
        strtime = time.strftime("%Y-%m-%d %H:%M:%S",crt)
        print(f"[{strtime}]",e)
        emit('flask_error',str(e))
        return
    emit('visData', new_check_data)
    with tempfile.TemporaryFile(mode="w+") as fp:
        python_write(map_data, fp)
        fp.seek(0)
        print("Send spdata")
        emit('genSpdata', fp.read())

if __name__ == "__main__":
    socketio.run(app,debug=False,host="0.0.0.0",port=5000)
