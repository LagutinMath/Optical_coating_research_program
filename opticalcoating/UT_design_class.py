from opticalcoating.design_class import Design
from json import load
# Инициализация

# ---------------------------------------------------
try:
    des = Design(name='Z36')
    print('OK')
except:
    print('Fail')

# ---------------------------------------------------
try:
    des = Design(d=[0., 100.], n_const=[0, 1.45])
    print('OK')
except:
    print('Fail')

# ---------------------------------------------------
try:
    des = Design(d=[0., 100.])
    print('Fail')
except:
    print('OK')

# ---------------------------------------------------
try:
    des = Design(d=[0., 100.])
    print('Fail')
except:
    print('OK')

# ---------------------------------------------------
with open('../Designs/Z36.json', 'r') as file:
    info = load(file)
try:
    des_HM44 = Design(info=info)
    print('OK')
except:
    print('Fail')

