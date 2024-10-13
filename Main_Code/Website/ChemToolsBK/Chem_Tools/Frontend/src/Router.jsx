import React, { useEffect } from 'react';
import { BrowserRouter as Router, Route, Routes, useLocation } from 'react-router-dom';
import App from './App.jsx';
import BalanceEquation from './pages/BalanceEquation';
import MolarityCalculator from './pages/MolarityCalculator';
import PHCalculator from './pages/PHCalculator';
import TitrationCurve from './pages/TitrationCurve';
import GibbsFreeEnergy from './pages/GibbsFreeEnergy';

const ScrollToTop = () => {
  const { pathname } = useLocation();

  useEffect(() => {
    window.scrollTo(0, 0);
  }, [pathname]);

  return null;
};

const AppRouter = () => {
  return (
    <Router>
      <ScrollToTop />
      <Routes>
        <Route path="/" element={<App />} />
        <Route path="/balance" element={<BalanceEquation />} />
        <Route path="/molarity" element={<MolarityCalculator />} />
        <Route path="/ph" element={<PHCalculator />} />
        <Route path="/titration" element={<TitrationCurve />} />
        <Route path="/gibbs" element={<GibbsFreeEnergy />} />
      </Routes>
    </Router>
  );
};

export default AppRouter;