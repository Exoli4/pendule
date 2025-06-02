from flask import Flask, request, jsonify

app = Flask(__name__)
CORS(app)

@app.route("/api/data", methods=["POST"])
def hello_world():
    return "<p>Hello, World!</p>"