def is_not_none(d: dict, k: str):
	if k in list(d.keys()):
		if d[k] is not None:
			return True
		else:
			return False
	else:
		return False
	