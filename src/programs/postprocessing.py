

@staticmethod
def poincare_hoff_rule(model):
    rule_text, rule_int = model.poincare_hoff_rule()
    if rule_int == 1:
        text = " is true!"
    else:
        text = " is not followed (" + str(rule_int) + ")"
    return "Poincare-Hoff rule " + rule_text + "= 1" + text
