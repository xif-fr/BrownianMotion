#include <SFML/Graphics.hpp>
#include <string>
#ifndef SFMLC01_WINDOW_HEIGHT
#define SFMLC01_WINDOW_HEIGHT SFMLC01_WINDOW_UNIT
#endif

namespace sf { namespace c01 {
	
	inline sf::Vector2f toWin (pt2_t p) {
		return sf::Vector2f( p.x*SFMLC01_WINDOW_UNIT, SFMLC01_WINDOW_HEIGHT-p.y*SFMLC01_WINDOW_UNIT );
	}
	
	inline void setCircleCR (sf::CircleShape& circle, pt2_t c, float r) {
		circle.setRadius(r*SFMLC01_WINDOW_UNIT);
		circle.setPosition( (c.x-r)*SFMLC01_WINDOW_UNIT, SFMLC01_WINDOW_HEIGHT-(c.y+r)*SFMLC01_WINDOW_UNIT );
	}
	inline sf::CircleShape buildCircleShapeCR (pt2_t c, float r) {
		sf::CircleShape circle; sf::c01::setCircleCR(circle, c, r); return circle;
	}
	
	inline void setRectShape (sf::RectangleShape& rect, pt2_t bottomleft, vec2_t diag) {
		rect.setSize(sf::Vector2f( diag.x*SFMLC01_WINDOW_UNIT, -diag.y*SFMLC01_WINDOW_UNIT ));
		rect.setPosition(sf::c01::toWin(bottomleft));
	}
	inline sf::RectangleShape buildRectShape (pt2_t bottomleft, vec2_t diag) {
		sf::RectangleShape rect; sf::c01::setRectShape(rect, bottomleft, diag); return rect;
	}
	
	inline sf::VertexArray buildLineStrip (const std::vector<pt2_t>& points, sf::Color color = sf::Color::Black) {
		sf::VertexArray lines (sf::LineStrip, points.size());
		for (size_t i = 0; i < points.size(); i++)
			lines[i] = sf::Vertex( sf::c01::toWin(points[i]), color );
		return lines;
	}
	inline sf::VertexArray buildLine (pt2_t a, pt2_t b, sf::Color color_a, sf::Color color_b) {
		sf::VertexArray lines (sf::Lines, 2);
		lines[0] = sf::Vertex( sf::c01::toWin(a), color_a );
		lines[1] = sf::Vertex( sf::c01::toWin(b), color_b );
		return lines;
	}
	inline sf::VertexArray buildLine (pt2_t a, pt2_t b, sf::Color color = sf::Color::Black) {
		return sf::c01::buildLine(a, b, color, color);
	}
	
	inline sf::ConvexShape buildParallelogram (pt2_t o, vec2_t va, vec2_t vb, sf::Color color) {
		sf::ConvexShape polygon;
		polygon.setPointCount(4);
		polygon.setPoint(0, sf::c01::toWin(o));
		polygon.setPoint(1, sf::c01::toWin(o+va));
		polygon.setPoint(2, sf::c01::toWin(o+va+vb));
		polygon.setPoint(3, sf::c01::toWin(o+vb));
		polygon.setFillColor(color);
		return polygon;
	}
	
	inline pt2_t fromWin (sf::Vector2f v) {
		return { .x = v.x/SFMLC01_WINDOW_UNIT, .y = 1.f-v.y/SFMLC01_WINDOW_UNIT };
	}
	
	inline sf::Text buildText (const sf::Font& font, pt2_t pos, const std::vector<std::wstring>& strs, sf::Color color = sf::Color::Black, uint8_t fontsize = 10) {
		sf::Text text;
		text.setFont(font);
		std::wstring s;
		for (auto& str : strs) {
			s += str;
			s += L"\n";
		}
		text.setString(s);
		text.setCharacterSize(fontsize);
		text.setFillColor(color);
		auto v = sf::c01::toWin(pos);
		text.setPosition((int)v.x, (int)v.y);
		return text;
	}

} }
