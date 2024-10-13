import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import '../App.css';

const MolarityCalculator = () => {
  const [solute, setSolute] = useState('');
  const [volume, setVolume] = useState('');
  const [molarity, setMolarity] = useState('');

  const handleSubmit = (e) => {
    e.preventDefault();
    // Funcionalidade para calcular a molaridade será implementada aqui
    setMolarity('A molaridade calculada aparecerá aqui');
  };

  return (
    <div className="ferramentas-quimica">
      <h1>Calculadora de Molaridade</h1>
      <p className="description">Ferramenta essencial para cálculos de concentração em soluções</p>
      
      <div className="conteudo">
        <div className="button-grid">
          <div className="category-column">
            <h2>Calculadora</h2>
            <form onSubmit={handleSubmit}>
              <div>
                <label htmlFor="solute">Digite a quantidade de soluto (moles):</label>
                <input
                  type="number"
                  id="solute"
                  value={solute}
                  onChange={(e) => setSolute(e.target.value)}
                  placeholder="ex., 1,5"
                  className="aba input-equation"
                />
              </div>
              <div>
                <label htmlFor="volume">Digite o volume da solução (litros):</label>
                <input
                  type="number"
                  id="volume"
                  value={volume}
                  onChange={(e) => setVolume(e.target.value)}
                  placeholder="ex., 1,0"
                  className="aba input-equation"
                />
              </div>
              <button type="submit" className="aba">Calcular Molaridade</button>
            </form>
            {molarity && (
              <div className="result">
                <h3>Molaridade:</h3>
                <p>{molarity}</p>
              </div>
            )}
          </div>
          <div className="category-column">
            <h2>Descrição do Problema</h2>
            <p>
              A molaridade é uma medida de concentração que expressa a quantidade de soluto
              em moles por litro de solução. Esta ferramenta ajuda você a calcular a molaridade
              de uma solução, fornecendo informações cruciais para experimentos e cálculos químicos.
            </p>
            <p>
              Para usar, insira a quantidade de soluto em moles e o volume da solução em litros
              nos campos apropriados. Clique em "Calcular Molaridade" e a ferramenta irá determinar
              a concentração molar da sua solução.
            </p>
          </div>
        </div>
      </div>
      <Link to="/" className="aba">Voltar para a Página Inicial</Link>
    </div>
  );
};

export default MolarityCalculator;
