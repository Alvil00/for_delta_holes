# coding:utf-8
import xml.etree.ElementTree as myxml


try:
	import matplotlib.font_manager as fm
	import matplotlib.pyplot as plt
	import matplotlib.patches as ptch
except ImportError:
	plt = None

try:
	import numpy as np
except ImportError:
	np = None
	print("Oops,routine can't find numpy! Try to connect with IT dep and install numpy module else you cant calculate the weakness ratio of your head!")

import random
import collections
import weakref
import math


if plt:
	GOST_FONT_FILE_NAME = "GOST.ttf"
	GOST_FONT = fm.FontProperties(fname=GOST_FONT_FILE_NAME)
	ISO_SYMBOLS_FONT_NAME = "ISOSymbols.ttf"
	ISO_SYMBOLS_FONT = fm.FontProperties(fname=ISO_SYMBOLS_FONT_NAME)


class AreaNode(object):
	def __init__(self, xleft: float, xright: float, ybot: float, ytop: float, hole_list, top_parent, deep):
		self._x = (xleft, xright)
		self._y = (ybot, ytop)
		self._list = []
		self._center = (xleft + xright) / 2, (ytop + ybot) / 2
		self._length = (xright - xleft), (ytop - ybot)
		self._areadict = None
		self._whole = False
		if top_parent is not None:
			self._top_parent = weakref.ref(top_parent)
		else:
			self._top_parent = None
		if hole_list:
			for i in hole_list:
				if isinstance(i, Hole):  # Hack with weak reference
					h = i
				else:
					h = i()
				if self.own(h):
					self._list.append(weakref.ref(h))
		if deep:
			self._areadict = dict()
			self._areadict[0] = AreaNode(self.center[0], xright, self.center[1], ytop, self._list, self.top_parent,
										 deep - 1)
			self._areadict[1] = AreaNode(xleft, self.center[0], self.center[1], ytop, self._list, self.top_parent,
										 deep - 1)
			self._areadict[2] = AreaNode(xleft, self.center[0], ybot, self.center[1], self._list, self.top_parent,
										 deep - 1)
			self._areadict[3] = AreaNode(self.center[0], xright, ybot, self.center[1], self._list, self.top_parent,
										 deep - 1)
		else:
			if self._top_parent is not None:
				self.top_parent._final_list.append(weakref.ref(self))

		if len(self._list) == 1:
			if self._areadict:
				if self._areadict[0]._whole and self._areadict[1]._whole and self._areadict[2]._whole and \
						self._areadict[3]._whole:
					self._whole = True
			else:
				a = self.point_in((self._x[0], self._y[0]))
				self._whole = a == self.point_in((self._x[0], self._y[1])) == self.point_in((self._x[1], self._y[0])) \
							  == self.point_in((self._x[1], self._y[1])) and a is not None

	@property
	def center(self):
		return self._center

	@property
	def top_parent(self):
		if self._top_parent is None:
			return None
		else:
			return self._top_parent()

	@property
	def length(self):
		return self._length

	def own(self, hole):
		_semilength = self.length[0] / 2, self.length[1] / 2
		_vec = (hole.pos[0] - self.center[0], hole.pos[1] - self.center[1])
		if abs(_vec[0]) <= _semilength[0] and abs(_vec[1]) <= _semilength[1]:
			return True
		if _vec[0] ** 2 + _vec[1] ** 2 <= (hole.diam / 2) ** 2:
			return True
		try:
			_tga = _vec[0] / _vec[1]
		except ZeroDivisionError:
			_tga = math.inf
		if _tga > 1.0 or _tga < -1.0 or _tga == math.inf:
			_vec = _vec[0] - math.copysign(_semilength[0], _vec[0]), _vec[1] * (1 - abs(_semilength[0] / _vec[0]))
			if abs(_vec[0]) <= hole.diam / 2:
				if self.center[1] - _semilength[1] <= hole.pos[1] <= self.center[1] + _semilength[1]:
					return True
			else:
				return False
		else:
			_vec = _vec[0] * (1 - abs(_semilength[1] / _vec[1])), _vec[1] - math.copysign(_semilength[1], _vec[1])
			if abs(_vec[1]) <= hole.diam / 2:
				if self.center[0] - _semilength[0] <= hole.pos[0] <= self.center[0] + _semilength[0]:
					return True
			else:
				return False
		if _vec[0] ** 2 + _vec[1] ** 2 <= (hole.diam / 2) ** 2:
			return True
		if _vec[0] > 0.0:
			_cx = self.center[0] + _semilength[0] - hole.pos[0]
		else:
			_cx = self.center[0] - _semilength[0] - hole.pos[0]
		if _vec[1] > 0.0:
			_cy = self.center[1] + _semilength[1] - hole.pos[1]
		else:
			_cy = self.center[1] - _semilength[1] - hole.pos[1]
		if _cx ** 2 + _cy ** 2 <= (hole.diam / 2) ** 2:
			return True
		return False

	# function for debugging of picture mode
	def point_in_square(self, pos):
		if self._areadict is None:
			return self
		else:
			if pos[0] >= self.center[0]:
				if pos[1] >= self.center[1]:
					return self._areadict[0].point_in_square(pos)
				else:
					return self._areadict[3].point_in_square(pos)
			else:
				if pos[1] >= self.center[1]:
					return self._areadict[1].point_in_square(pos)
				else:
					return self._areadict[2].point_in_square(pos)

	def point_in(self, pos):
		if not self._list:
			return None
		if self._areadict is None:
			if self._whole:
				return self._list[-1]()
			for i in self._list:
				if (pos[0] - i().pos[0]) ** 2 + (pos[1] - i().pos[1]) ** 2 <= (i().diam / 2) ** 2:
					return i()
		else:
			if self._whole:
				return self._list[-1]()
			if pos[0] >= self._center[0]:
				if pos[1] >= self._center[1]:
					return self._areadict[0].point_in(pos)
				else:
					return self._areadict[3].point_in(pos)
			else:
				if pos[1] >= self._center[1]:
					return self._areadict[1].point_in(pos)
				else:
					return self._areadict[2].point_in(pos)
		return None

	def has_intersection_with_rectangle(self, xleft: float, xright: float, ybot: float, ytop: float) -> bool:
		if xright < self._x[0] or xleft > self._x[1] or ytop < self._y[0] or ybot > self._y[1]:
			return False
		if not self._list:
			return False
		if self._whole:
			return True
		args = xleft, xright, ybot, ytop
		if self._areadict is None:
			temp_areanode = AreaNode(*args, None, None, 0)
			for hole in self._list:
				if temp_areanode.own(hole()):
					return True
			return False
		else:
			result = False
			for i in range(4):
				result = result or self._areadict[i].has_intersection_with_rectangle(*args)
				if result:
					return True


class AreaTree(AreaNode):
	def __init__(self, xleft, xright, ybot, ytop, hole_list, deep=7):
		self._final_list = []
		super().__init__(xleft, xright, ybot, ytop, hole_list, self, deep)

	def __getitem__(self, pos):
		return super().point_in(pos)

	def walk_final(self):
		for i in self._final_list:
			yield i()


class Plita(object):
	WEDGE_THETA_ANGLE_STYLE = (
		(90.0, 180.0),
		(-90.0, 90.0),
		(0.0, 90.0),
		(180.0, 270.0),
		(0.0, 180.0),
		(90.0, 270.0),
		(270.0, 360.0),
		(180.0, 360.0)
	)

	def __init__(self, diam, thickness, thinning, sigma, pressure):
		self.d = diam
		self.s = thickness
		self._types_of_holes = dict()
		if isinstance(thinning, collections.Iterable):
			if len(thinning) == 3:
				self.thinning = thinning
			else:
				raise ValueError("Wrong number of thinning args")
		elif type(thinning) in (float, int):
			self.thinning = thinning
		else:
			raise ValueError("Wrong type of thinning")

		if isinstance(sigma, collections.Iterable):
			if len(sigma) == 2:
				self.sigma = int(min(max(sigma) / 2.6, min(sigma) / 1.5))
			else:
				raise ValueError("Wrong number of sigma args")
		elif type(sigma) in (float, int):
			self.sigma = sigma
		else:
			raise ValueError("Wrong type of sigma")

		self.p = pressure
		self.holes = None
		self._fi = 0.0
		self._fi_ang = 0.0
		self.s_i = 0.0
		self.s_r = 0.0
		self._areatree = None

	def _append_holes(self, holes):
		self.holes = holes

	def append_holes_by_types(self, holes_diams: dict, holes_positions: dict):
		for i in holes_diams.keys():
			Hole.create_type(i, holes_diams.get(i))
			if not self._types_of_holes.get(i):
				self._types_of_holes[i] = Hole.__dict__[i]
		hole_list = []
		for i in holes_diams.keys():
			pos = holes_positions.get(i)
			if pos:
				for cpos in pos:
					hole_list.append(Hole.__dict__[i](*cpos))
		self._append_holes(hole_list)

	def create_doc(self, filename=None):
		if filename is None:
			filename = "output.xml"
		main_elem = myxml.Element("inputData")
		myxml.SubElement(main_elem, "d", {"value": str(self.d - 1).replace(".", ",")})
		myxml.SubElement(main_elem, "d2", {"value": str(self.d).replace(".", ",")})
		myxml.SubElement(main_elem, "s1", {"value": str(self.s).replace(".", ",")})
		myxml.SubElement(main_elem, "s2", {"value": str(self.s).replace(".", ",")})
		if type(self.thinning) in (list, tuple):
			a = self.thinning
		else:
			a = [self.thinning, 0, 0]
		myxml.SubElement(main_elem, "c11", {"value": str(a[0]).replace(".", ",")})
		myxml.SubElement(main_elem, "c12", {"value": str(a[1]).replace(".", ",")})
		myxml.SubElement(main_elem, "c2", {"value": str(a[2]).replace(".", ",")})
		myxml.SubElement(main_elem, "sigma", {"value": str(self.sigma).replace(".", ",")})
		myxml.SubElement(main_elem, "P", {"value": str(self.p).replace(".", ",")})
		if self.holes is not None:
			ta = myxml.SubElement(main_elem, "tableApertures")
			for i in self.holes:
				current_hole = myxml.SubElement(ta, "string")
				current_hole.text = "{diam} {posx} {posy}".format(diam=i.diam, posx=i.pos[0], posy=i.pos[1]).replace(
					".", ",")

		myxml.ElementTree(element=main_elem).write(filename, encoding="utf-8", xml_declaration=True)

	# debugging function for pick test
	def _pick_test(self, event):
		self._areatree.point_in_square((event.xdata, event.ydata))._rec.set_color('g')
		print(self._areatree.point_in_square((event.xdata, event.ydata))._center)

	def save_pic(self, filename, weak_line=True, to_print=True, dimensions=True, wedges=True, holenames=True,
				 dimension_random_position=True, debug=False) -> None:
		"""Function to save picture of the result"""

		def box_intersection(box1: tuple, box2: tuple) -> bool:
			if ((box1[0] <= box2[0] <= box1[1] or box1[0] <= box2[1] <= box1[1]) or \
				(box2[0] < box1[0] < box2[1] or box2[0] < box1[1] < box2[1])) and \
					((box1[2] <= box2[2] <= box1[3] or box1[2] <= box2[2] <= box1[3]) or \
					 (box2[2] < box1[2] < box2[3] or box2[2] < box1[2] < box2[3])):
				return True

		def out_of_plita_range(box:tuple) -> bool:
			nonlocal radius
			max_abs_x_coord = max(math.fabs(box[0]), math.fabs(box[1]))
			max_abs_y_coord = max(math.fabs(box[2]), math.fabs(box[3]))
			if max_abs_x_coord ** 2 + max_abs_y_coord ** 2 > lim ** 2:
				return True
			else:
				return False

		def print_iso_dimension_annotation(annotate_hole: Hole, angle: float = 1.5 * math.pi / 4,
										   size: float = 0.15) -> None:
			"""Function to annotate the dimension of the hole"""

			nonlocal ax
			nonlocal lim
			nonlocal _fh
			nonlocal _text_boxes
			nonlocal weak_line
			dimension_text = "Ø{:<5.1f}".format(annotate_hole.diam).replace(".", ",")
			# Удаляем 2 последних символа (если это цел. диаметр), чтобы писать размеры близкие к целочисл. без запятой
			if dimension_text.strip()[-1] == '0':
				dimension_text = dimension_text.strip()[:-2]
			# Блок вычисления границ аннотации и текстбокса
			factor = math.cos(angle), math.sin(angle)
			pan_length = 0.133 * lim
			an_pos = annotate_hole.pos[0] + annotate_hole.diam / 2 * factor[0], \
					 annotate_hole.pos[1] + annotate_hole.diam / 2 * factor[1]
			an_text_pos = an_pos[0] + size * lim * factor[0], \
						  an_pos[1] + size * lim * factor[1]
			text_box = an_text_pos[0] - lim * 0.01, an_text_pos[0] + pan_length, \
					   an_text_pos[1] - lim * 0.015, an_text_pos[1] + lim * 0.07
			# Проверка не пересекает ли текстбокс отверстие
			if self._areatree.has_intersection_with_rectangle(*text_box):
				return None
			# Проверка не пересекает ли данный текстбокс другие текстбоксы
			if _text_boxes:
				for tb in _text_boxes:
					if box_intersection(text_box, tb):
						return None
			# Проверка по контрольным точкам в стрелке
			num_of_control_points = 5
			control_point_semi_side = lim * 0.002
			cp_box_list = list()
			for i in range(1, num_of_control_points + 1):
				cp_center_pos = (i * (an_text_pos[0] - an_pos[0]) / (num_of_control_points + 1) + an_pos[0],
								 i * (an_text_pos[1] - an_pos[1]) / (num_of_control_points + 1) + an_pos[1])
				cp_box = (cp_center_pos[0] - control_point_semi_side, cp_center_pos[0] + control_point_semi_side,
						  cp_center_pos[1] - control_point_semi_side, cp_center_pos[1] + control_point_semi_side)
				if self._areatree.has_intersection_with_rectangle(*cp_box):
					return None
				for tb in _text_boxes:
					if box_intersection(cp_box, tb):
						return None
				cp_box_list.append(cp_box)
			# Проверка не лежит ли текстбокс за пределами плиты
			if out_of_plita_range(text_box):
				return None
			# Проверка не попадает ли текстбокс на горизонтальную ось
			if text_box[3] > 0.0 > text_box[2]:
				return None
			# Проверка не попадает ли текстбокс на линию минимального сечения
			if weak_line:
				k = math.tan(self._fi_ang)
				if (math.isfinite(k) or math.fabs(k) <= 60.0) and \
						(not math.copysign(1.0, text_box[0] * k - text_box[2]) +
							 math.copysign(1.0, text_box[1] * k - text_box[3]) or
						 not math.copysign(1.0, text_box[0] * k - text_box[3]) +
							 math.copysign(1.0, text_box[1] * k - text_box[2])):
					return None
			# По итогу всех проверок вставляем текстбоксы в общий словарь текстбоксов
			_text_boxes.append(text_box)
			_text_boxes.extend(cp_box_list)
			# Если размер попадает на вретикальную ось симметрии плиты окрашиваем участок который он проходит белым
			if text_box[1] > 0.0 > text_box[0]:
				ax.add_artist(plt.Line2D((0.0, 0.0),
										 (an_text_pos[1] - lim * 0.004, an_text_pos[1] + lim * 0.057),
										 color='white', linewidth=0.35))
			# Далее рисуем аннотацию без линии стрелки
			plt.annotate(dimension_text, xy=an_pos, xytext=an_text_pos,
						 arrowprops=dict(arrowstyle='-|>', facecolor='black', linewidth=0.25,
										 relpos=(0.0, 0.0), shrinkA=-12.0, shrinkB=0.0),
						 fontsize=8, fontproperties=ISO_SYMBOLS_FONT, fontweight="light", va="bottom")
			# Рисуем полку
			ax.add_artist(plt.Line2D((an_text_pos[0], an_text_pos[0] + pan_length),
									 (an_text_pos[1] - lim * 0.002, an_text_pos[1] - lim * 0.002),
									 color='black', linewidth=0.25))
			# Рисуем линию стрелки
			ax.add_artist(plt.Line2D((an_pos[0], an_text_pos[0]),
									 (an_pos[1], an_text_pos[1] - lim * 0.002),
									 color='black', linewidth=0.25))
			# Ставим флаг что размер проставлен
			_fh[type(annotate_hole)] = True

		def print_iso_name_annotation(annotate_hole, angle: float, size:float, inner_sizes: bool = True,
									  debug: bool = False):
			nonlocal ax
			nonlocal lim
			nonlocal _text_boxes
			nonlocal max_hole_type
			if inner_sizes and isinstance(hole, max_hole_type) and max_hole_type._STD_DIAM >= 0.175 * lim:
				_sm = 0.56
			else:
				_sm = 1.0
			if annotate_hole.name is None:
				name = type(annotate_hole).__name__
			else:
				name = annotate_hole.name
			factor = math.cos(angle), math.sin(angle)
			text_point = annotate_hole.pos[0] + (_sm * annotate_hole.diam / 2 + size * lim) * factor[0], \
						 annotate_hole.pos[1] + (_sm * annotate_hole.diam / 2 + size * lim) * factor[1]
			if factor[0] >= 0.0:
				tb_x = text_point[0] - 0.025 * lim, text_point[0] + 0.0195 * lim * len(name)
				_m = 1.0
			else:
				tb_x = text_point[0] - 0.0195 * lim * len(name) - 0.022 * lim, text_point[0]
				_m = 0.1
			if factor[1] >= 0.0:
				tb_y = text_point[1] - 0.0125 * lim, text_point[1] + 0.0275 * lim
			else:
				tb_y = text_point[1] - 0.03 * lim, text_point[1] + 0.0125 * lim
			text_box = *tb_x, *tb_y
			if not (inner_sizes and isinstance(hole, max_hole_type)):
				if self._areatree.has_intersection_with_rectangle(text_box[0]-0.075 * lim * _m, *text_box[1:]):
					return None
			for tb in _text_boxes:
				if box_intersection(text_box, tb):
					return None
			if out_of_plita_range(text_box):
				return None
			if text_box[3] > 0.0 > text_box[2]:
				return None
			if debug:
				ax.add_artist(plt.Rectangle((text_box[0], text_box[2]),
											(text_box[1] - text_box[0]), (text_box[3] - text_box[2]),
											color='red', alpha=0.3))
			if text_box[1] > 0.0 > text_box[0]:
				ax.add_artist(plt.Line2D((0.0, 0.0),
										 (text_box[2] - lim * 0.011, text_box[3] + lim * 0.01),
										 color='white', linewidth=0.35))
			_text_boxes.append(text_box)
			ax.add_artist(plt.Text(text_box[0],text_box[2], name,
								   fontsize=6, fontproperties=GOST_FONT, fontweight="light"))
			return True

		if plt:
			radius = self.d / 2
			lim = radius
			fig, ax = plt.subplots()
			if to_print:
				# Тип штрихпунктирной линии (дает красивые результаты только для центральных линий)
				lstyle = (-10, (20, 4, 4, 4))
				# Для отвестия с максимальным диаметром предполагается пустая штриховка
				max_hole_type = max(list(self._types_of_holes.items()), key=lambda a: a[1]._STD_DIAM)[1]
				c1 = plt.Circle((0, 0), radius, color='black', fill=False, linewidth=0.5)
				ax.add_artist(c1)
				# Горизонтальая ось симметрии плиты
				ax.add_artist(
					plt.Line2D((-1.05 * lim, 1.05 * lim), (0.0, 0.0), linestyle=lstyle, linewidth=0.3, color='black'))
				# Вертикальная ось симметрии плиты
				ax.add_artist(
					plt.Line2D((0.0, 0.0), (-1.05 * lim, 1.05 * lim), linestyle=lstyle, linewidth=0.3, color='black'))
				# Штриховка - назначение типов для отверстий
				wedge_iter = iter(self.WEDGE_THETA_ANGLE_STYLE)
				_hd = dict()  # Словарь содержит типы обозначения-штриховки отверстий
				_fh = dict()  # Словарь флагиов,которые отвечают за то обозначен ли уже размер отверстия для даного типа
				_text_boxes = list()  # Список текстбоксов которые находятся на холсте
				for hole_type in self._types_of_holes.keys():
					if max_hole_type is not self._types_of_holes[hole_type]:
						_hd[self._types_of_holes[hole_type]] = next(wedge_iter)
					else:
						_hd[self._types_of_holes[hole_type]] = None
					_fh[self._types_of_holes[hole_type]] = False
				standard_annotation_length = 0.15  # Длинна стрелки в размере отверстия
				annotation_lengthes = (
					standard_annotation_length * 0.85, standard_annotation_length, standard_annotation_length * 1.1)
				annotation_name_lengthes = (0.021, 0.03, 0.04, 0.05, 0.064, 0.07)
				s = 4  # Переменная влияющая на количество уговых точек
				ang_dispersion = lambda x, s=s: \
					[(x / (2 * s) - 0.01) * math.pi, x / (2 * s) * math.pi, (x / (2 * s) + 0.01) * math.pi]
				_ad = [ang_dispersion(x) for x in range(-s, 2 * s) if x % s != 0]
				s = 24
				_fn = [ang_dispersion(x) for x in range(-2 * s, 2 * s) if x != -2 * s]
				if dimension_random_position:
					_ad = random.sample(_ad, len(_ad))
				annotation_angles = list()
				annotation_name_angles = list()
				for te_ in _fn:
					annotation_name_angles.extend(te_)
				for te_ in _ad:
					annotation_angles.extend(te_)
			else:
				c1 = plt.Circle((0, 0), radius, color='r', fill=False)
				ax.add_artist(c1)
			for hole in self.holes:
				if to_print:
					hr = hole.diam / 2  # Радиус отверстия
					hlr = hr * 1.075 + radius * 0.01  # Половина длинны осей симметрии отверстия
					# Если радиус отверстия состовляет менее 8 % от радиуса плиты - тип линий оси будует рисоваться
					# сплошной линией иначе штрихпунктирной
					if hr >= radius * 0.08:
						hls = (-30, lstyle[1])
					else:
						hls = None
					# Рисуем отверстие
					ax.add_artist(plt.Circle(hole.pos, hr, color='black', fill=False, linewidth=0.5))
					# Вычисляем границы для бокса линий симметрии отверстия
					hlim = hole.pos[0] - hlr, hole.pos[0] + hlr, hole.pos[1] - hlr, hole.pos[1] + hlr
					# Горизонтальная ось симметрии для отверстий
					ax.add_artist(
						plt.Line2D((hlim[0], hlim[1]), (hole.pos[1], hole.pos[1]), linestyle=hls, linewidth=0.3,
								   color='black'))
					# Вертикальная ось симметрии для отверстий
					ax.add_artist(
						plt.Line2D((hole.pos[0], hole.pos[0]), (hlim[2], hlim[3]), linestyle=hls, linewidth=0.3,
								   color='black'))
					# Добавляем штриховку к отверстию если в словаре есть штриховка для соотв. типа
					if _hd[type(hole)] and wedges:
						ax.add_artist(ptch.Wedge(hole.pos, hr, *_hd[type(hole)], fill=True, facecolor='black'))
					# Пытаемся добавить размер отверстия если это необходимо
					if not _fh[type(hole)] and dimensions:
						for annotation_length in annotation_lengthes:
							for annotation_angle in annotation_angles:
								print_iso_dimension_annotation(hole, annotation_angle, annotation_length)
								if _fh[type(hole)]: break
							if _fh[type(hole)]: break
				else:
					ax.add_artist(plt.Circle(hole.pos, hole.diam / 2, color='g', fill=False))
			if to_print and holenames:
				for hole in self.holes:
					tanno = None
					for annotation_length in annotation_name_lengthes:
						for annotation_angle in annotation_name_angles:
							tanno = print_iso_name_annotation(hole, annotation_angle, annotation_length)
							if tanno is not None:
								break
						if tanno is not None:
							break

			ax.set_xlim((-1.1 * lim, 1.1 * lim))
			ax.set_ylim((-1.1 * lim, 1.1 * lim))
			plt.axis('scaled')
			if weak_line:
				weak_line_artist = plt.Line2D(self._nline[0], self._nline[1], linewidth=1)
				ax.add_artist(weak_line_artist)
			if debug:
				for num, area in enumerate(self._areatree.walk_final()):
					area._rec = plt.Rectangle((area._x[0], area._y[0]), area.length[0], area.length[1])
					if area._list:
						area._rec.set_color('b')
						if area._whole:
							area._rec.set_color('y')
					else:
						area._rec.set_color('r')
					area._rec.set_alpha(0.2)
					ax.add_artist(area._rec)
					fig.canvas.mpl_connect('button_press_event', self._pick_test)
			plt.axis([-1.1 * lim, 1.1 * lim, -1.1 * lim, 1.1 * lim])
			if debug:
				plt.show()
			else:
				plt.savefig(filename, dpi=600)
		else:
			print("no mplt")

	def factorise(self, deep=7):
		self._areatree = AreaTree(-self.d / 2, self.d / 2, -self.d / 2, self.d / 2, self.holes, deep)

	def calc_fi_geom(self, start_angle=0.0, final_angle=180.0, dim=1, dimfi=1):
		if np:
			num_of_steps_linear = int(self.d) * dim
			total_angle = final_angle - start_angle
			num_of_steps_angle = math.ceil(dimfi * total_angle) + 1
			fi_max = 1.0
			d = 0.0
			d_max = 0.0
			fi_max_angle = 0.0
			int_d = np.linspace(-self.d / 2, self.d / 2, num_of_steps_linear, dtype=np.float64)
			step_ = self.d / num_of_steps_linear
			ind_f = np.linspace(math.radians(start_angle), math.radians(final_angle), num_of_steps_angle,
								dtype=np.float64)
			for f in ind_f:
				coords1 = np.array([(np.cos(f) * int_d), (np.sin(f) * int_d)])
				coords2 = coords1.swapaxes(0, 1)
				for x, y in coords2:
					if self._areatree[(x, y)] is not None:
						d += step_
				k = (d / self.d)
				fi_ = 1 / (1 + k + k * k)
				if fi_ < fi_max:
					self._nline = coords1
					fi_max = fi_
					fi_max_angle = f
					d_max = d
				d = 0.0
			self._fi = fi_max
			self._fi_ang = fi_max_angle
			return fi_max, math.degrees(fi_max_angle), d_max
	
	def calc_fi_analytic(self, start_angle=0.0, final_angle=180.0, dim=1, dimfi=1):
		if np:
			total_angle = final_angle - start_angle
			num_of_steps_angle = math.ceil(dimfi * total_angle) + 1
			fi_max = 1.0
			d = 0.0
			d_max = 0.0
			fi_max_angle = 0.0
			ind_f = np.linspace(math.radians(start_angle), math.radians(final_angle), num_of_steps_angle,
								dtype=np.float64)
			for f in ind_f:
				d = 0.0
				for hole in self.holes:
					#ay**2+by+c = 0
					a = math.tan(f) ** 2 + 1.0
					b = - 2.0 * (- hole.pos[1] * math.tan(f) + hole.pos[0])
					c = hole.pos[0] ** 2 + hole.pos[1] ** 2 - (hole.diam /2) ** 2
					discr = b ** 2 - 4 * a * c
					if discr > 0.0:
						discr = math.sqrt(discr)
						y = (- b + discr ) / (2 * a), (- b - discr ) / (2 * a)
						x = y[0] * math.tan(f), y[1] * math.tan(f)
						d += math.sqrt((y[0] - y[1]) ** 2 + (x[0] - x[1]) ** 2)
					else:
						continue
					k = (d / self.d)
					fi_ = 1 / (1 + k + k * k)
				if fi_ < fi_max:
					fi_max = fi_
					fi_max_angle = f
					d_max = d
			int_d = np.linspace(-self.d / 2, self.d / 2, 2, dtype=np.float64)
			self._nline = np.array([(np.cos(fi_max_angle) * int_d), (np.sin(fi_max_angle) * int_d)])
			self._fi = fi_max
			self._fi_ang = fi_max_angle
			return fi_max, math.degrees(fi_max_angle), d_max
	
	def calc_s(self):
		self.s_i = 0.45 * self.d * math.sqrt(self.p / (self._fi * self.sigma))
		if type(self.thinning) in (tuple, list):
			self.s_r = self.s_i + math.fsum(self.thinning)
		else:
			self.s_r = self.s_i + self.thinning
		return self.s, self.s_i, self.s_r


class Hole:
	@classmethod
	def create_type(cls, typename, diam):
		def _create(self, posx, posy, hole_name=None):
			return self.__class__.__mro__[1].__init__(self, (posx, posy), self._STD_DIAM, hole_name)

		cls.__class__.__setattr__(cls, typename, type(typename, (cls,), {'_STD_DIAM': diam, '__init__': _create}))

	def __init__(self, pos: tuple, diam: float, name: str = None):
		if not isinstance(pos, tuple):
			raise type("TupleError", (ValueError,),
					   {"__init__": lambda self, a: print("wrong coords type: {}\n".format(a))})(type(pos))
		if len(pos) != 2:
			raise type("CoordError", (ValueError,),
					   {"__init__": lambda self, a: print("wrong coords num: {}".format(a))})(pos)
		self.pos = pos
		self.diam = diam
		self.name = name

	def __str__(self):
		r = "<string> {diam} {posx} {posy} </string>\n"
		return r.format(diam=self.diam, posx=self.pos[0], posy=self.pos[1])


# Основная функция расчета fi_calc_macro
# args:
# p - экземпляр класса Plita для которго нужен расчет;
# filename - имя файла графического представления плиты(будет создан);
# fn - степень факторизации (оптимально равна 7,можно увеличить при наличии большого количества маленьких отверстий);
# acc_lin - линейная точность (Для выполнению расчета достаточно точности вплоть до 1.0,для кросс-верификации
# с ПК DELTA можно взять точность вплоть до 0.1).


def fi_calc_macro(p: Plita, filename: str = "fig.png", fn: int = 7, acc_lin: float = 1.0, sp_kwargs: dict = None, method=Plita.calc_fi_analytic) \
		-> None:
	if np:
		if acc_lin <= 0.0:
			raise ValueError("Linear accuracy must be greater than zero.")
		if sp_kwargs is None:
			sp_kwargs = dict()
		if sp_kwargs.get("debug"):
			sp_kwargs.pop("debug")
		acc_lin_ = math.ceil(1 / acc_lin)
		print("factorisation...")
		p.factorise(fn)
		print("calculating fi..")
		preliminary_angle = method(p)[1]
		a = (math.floor(preliminary_angle) - 0.5, math.ceil(preliminary_angle) + 0.5)
		print("fi:{res[0]:<7.5f}\tangle_deg:{res[1]:<6.2f}\tlength:{res[2]:<8.1f}".format(
			res=method(p, a[0], a[1], 2 * acc_lin_, 16)))
		print("s_assumption:{res[0]:<7.3f}\ns_calculated:{res[1]:<7.3f}\ns_total:{res[2]:<7.3f}".format(res=p.calc_s()))
		p.save_pic(filename, **sp_kwargs)


def example():
	holes_positions = dict()
	p = Plita(diam=2190, thickness=115, thinning=(0.35, 0, 0), sigma=124, pressure=0.6)
	hole_type_diameters = {
		'KG': (22.0 * 131.0 + 56.0 * 122.0 + 37.0 * 127.0) / (22.0 + 56.0 + 37.0),
		'PR': (22.0 * 90.0 + 81.0 * 81.0 + 12.0 * 88.0) / (22.0 + 81.0 + 12.0),
		'TS': (22.0 * 58.0 + 81.0 * 49.0 + 12.0 * 56.0) / (22.0 + 81.0 + 12.0),
		'T': (22.0 * 46.0 + 81.0 * 37.0 + 12.0 * 44.0) / (22.0 + 81.0 + 12.0),
		'EG': (22.0 * 58.0 + 81.0 * 49.0 + 12.0 * 56.0) / (22.0 + 81.0 + 12.0),
		'GK': (15.0 * 608.0 + 20.0 * 570.0 + 80.0 * 563.0) / (15.0 + 20.0 + 80.0)
	}
	# List of GK holes
	holes_positions["GK"] = [
		(402.94, 575.45,   "ГК-1"),  # GK1
		(-402.94, -575.45, "ГК-2"),  # GK2
		(575.45, -402.94,  "ГК-3"),  # GK3
		(-575.45, 402.94,  "ГК-4")  # GK4
	]
	# List of KG holes
	holes_positions["KG"] = [
		(-151.5, 29.16,     "ЦКГ-1"),  # CKG1
		(151.5, -29.16,     "ЦКГ-2"),  # CKG2
		(-303.0, -174.94,   "СКГ-1"),  # SKG1
		(-151.5, 320.72,    "СКГ-2"),  # SKG2
		(303.0, 174.94,     "СКГ-3"),  # SKG3
		(151.5, -320.72,    "СКГ-4"),  # SKG4
		(-555.5, -87.47,    "ПКГ-1"),  # PKG1
		(0.0, 524.81,       "ПКГ-2"),  # PKG2
		(555.5, 87.47,      "ПКГ-3"),  # PKG3
		(0.0, -524.81,      "ПКГ-4")  # PKG4
	]
	# List of TS holes
	holes_positions["TS"] = [
		(-333.5, 132.47,    "ТС1"),  # TS1
		(62.03, 327.37,     "ТС2"),  # TS2
		(89.47, 109.97,     "ТС3"),  # TS3
		(151.5, -482.34,    "ТС4"),  # TS4
		(-101.0, -219.94,   "ТС5")  # TS5
	]
	# List of PR holes
	holes_positions["PR"] = [
		(30.1, 230.0,   "Пр1"),  # PR1
		(-29.8, -130.1, "Пр2")  # PR2
	]
	# List of EG holes
	holes_positions["EG"] = [
		(-252.5, 612.28, "ЭГ1"),  # EG1
		(252.5, -612.28, "ЭГ2")   # EG2
	]
	# List of T holes
	holes_positions["T"] = [
		(-521.0, 744.0, "Твх.1"),  # Tvh1
		(744.0, 521.0,  "Твх.2"),  # Tvh2
		(521.0, -744.0, "Твх.3"),  # Tvh3
		(-744.0, -521.0, "Твх.4"),  # Tvh4
		(-270.0, 189.0, "Твых.1"),  # Tvyh1
		(165.0, 285.0, "Твых.2"),  # Tvyh2
		(270.0, -189.0, "Твых.3"),  # Tvyh3
		(-165.0, -285.0, "Твых.4")  # Tvyh4
	]
	p.append_holes_by_types(hole_type_diameters, holes_positions)
	p.create_doc("list.xml")
	fi_calc_macro(p)


if __name__ == "__main__":
	example()
